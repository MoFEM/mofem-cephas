/** \file ElementsOnEntities.cpp

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

#include <Includes.hpp>
#include <version.h>
#include <definitions.h>
#include <Common.hpp>

#include <h1_hdiv_hcurl_l2.h>
#include <fem_tools.h>

#include <UnknownInterface.hpp>

#include <MaterialBlocks.hpp>
#include <CubitBCData.hpp>
#include <TagMultiIndices.hpp>
#include <CoordSysMultiIndices.hpp>
#include <FieldMultiIndices.hpp>
#include <EntsMultiIndices.hpp>
#include <DofsMultiIndices.hpp>
#include <FEMMultiIndices.hpp>
#include <ProblemsMultiIndices.hpp>
#include <AdjacencyMultiIndices.hpp>
#include <BCMultiIndices.hpp>
#include <CoreDataStructures.hpp>
#include <SeriesMultiIndices.hpp>

#include <LoopMethods.hpp>
#include <FieldInterface.hpp>
#include <MeshRefinment.hpp>
#include <PrismInterface.hpp>
#include <SeriesRecorder.hpp>
#include <Core.hpp>

#include <FTensor.hpp>
#include <DataStructures.hpp>
#include <DataOperators.hpp>
#include <ElementsOnEntities.hpp>

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

static MoABErrorCode rval;

PetscErrorCode ForcesAndSurcesCore::getNumberOfNodes(int &num_nodes) const {
  PetscFunctionBegin;

  EntityHandle ent = numeredEntFiniteElementPtr->get_ent();
  switch(mField.get_moab().type_from_handle(ent)) {
    case MBVERTEX:
    num_nodes = 1;
    break;
    case MBEDGE:
    num_nodes = 2;
    break;
    case MBTRI:
    num_nodes = 3;
    break;
    case MBTET:
    num_nodes = 4;
    break;
    case MBPRISM:
    num_nodes = 6;
    break;
    default:
    SETERRQ(PETSC_COMM_SELF,MOFEM_NOT_IMPLEMENTED,"not implemented");
  }

  PetscFunctionReturn(0);
}

// ** Sense **

PetscErrorCode ForcesAndSurcesCore::getSense(EntityType type,boost::ptr_vector<DataForcesAndSurcesCore::EntData> &data) const {
  PetscFunctionBegin;
  try {
    SideNumber_multiIndex& side_table = const_cast<SideNumber_multiIndex&>(numeredEntFiniteElementPtr->get_side_number_table());
    if(data.size() < side_table.get<2>().count(type)) {
      // prims has 9 edges, some of edges for "flat" prism are not active
      SETERRQ2(PETSC_COMM_SELF,MOFEM_DATA_INCONSISTENCY,"data inconsistency %u < %u",data.size(),side_table.get<2>().count(type));
    }
    SideNumber_multiIndex::nth_index<2>::type::iterator siit = side_table.get<2>().lower_bound(type);
    SideNumber_multiIndex::nth_index<2>::type::iterator hi_siit = side_table.get<2>().upper_bound(type);
    for(;siit!=hi_siit;siit++) {
      data[siit->get()->side_number].getSense() = siit->get()->sense;
      if(siit->get()->brother_side_number!=-1) {
        if(data.size() < (unsigned)siit->get()->brother_side_number) {
          SETERRQ(PETSC_COMM_SELF,MOFEM_DATA_INCONSISTENCY,"data inconsistency");
        }
        data[siit->get()->brother_side_number].getSense() = siit->get()->sense;
      }
    }
  } catch (std::exception& ex) {
    std::ostringstream ss;
    ss << "thorw in method: " << ex.what() << " at line " << __LINE__ << " in file " << __FILE__;
    SETERRQ(PETSC_COMM_SELF,MOFEM_STD_EXCEPTION_THROW,ss.str().c_str());
  }

  PetscFunctionReturn(0);
}

PetscErrorCode ForcesAndSurcesCore::getEdgesSense(DataForcesAndSurcesCore &data) const {
  PetscErrorCode ierr;
  PetscFunctionBegin;
  ierr = getSense(MBEDGE,data.dataOnEntities[MBEDGE]); CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

PetscErrorCode ForcesAndSurcesCore::getTrisSense(DataForcesAndSurcesCore &data) const {
  PetscErrorCode ierr;
  PetscFunctionBegin;
  ierr = getSense(MBTRI,data.dataOnEntities[MBTRI]); CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

PetscErrorCode ForcesAndSurcesCore::getQuadSense(DataForcesAndSurcesCore &data) const {
  PetscErrorCode ierr;
  PetscFunctionBegin;
  ierr = getSense(MBQUAD,data.dataOnEntities[MBQUAD]); CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

// ** Order **

template<typename DOFMULTIINDEX>
static int get_max_order(
  const DOFMULTIINDEX &dof_multi_index
) {
  int max_order = 0;
  typename DOFMULTIINDEX::iterator dit,hi_dit;
  dit = dof_multi_index.begin();
  hi_dit = dof_multi_index.end();
  for(;dit!=hi_dit;dit++) {
    int dit_max_order = (*dit)->get_max_order();
    max_order = (max_order>dit_max_order) ? max_order : dit_max_order;
  }
  return max_order;
}

int ForcesAndSurcesCore::getMaxDataOrder() const {
  return get_max_order(numeredEntFiniteElementPtr->get_data_dofs());
}

int ForcesAndSurcesCore::getMaxRowOrder() const {
  return get_max_order(numeredEntFiniteElementPtr->get_rows_dofs());
}

int ForcesAndSurcesCore::getMaxColOrder() const {
  return get_max_order(numeredEntFiniteElementPtr->get_cols_dofs());
}

PetscErrorCode ForcesAndSurcesCore::getDataOrder(const EntityType type,const FieldSpace space,boost::ptr_vector<DataForcesAndSurcesCore::EntData> &data) const {
  PetscFunctionBegin;
  try {
    SideNumber_multiIndex& side_table = const_cast<SideNumber_multiIndex&>(numeredEntFiniteElementPtr->get_side_number_table());
    if(data.size() < side_table.get<2>().count(type)) {
      // prims has 9 edges, some of edges for "flat" prism are not active
      SETERRQ2(PETSC_COMM_SELF,MOFEM_DATA_INCONSISTENCY,
        "data inconsistency %d < %d",
        data.size(),side_table.get<2>().count(type)
      );
    }
    for(unsigned int side = 0;side<data.size();side++) {
      data[side].getDataOrder() = 0;
    }
    FEDofEntity_multiIndex::index<Composite_EntType_and_Space_mi_tag>::type &data_dofs =
    const_cast<FEDofEntity_multiIndex::index<Composite_EntType_and_Space_mi_tag>::type&>(
      numeredEntFiniteElementPtr->get_data_dofs().get<Composite_EntType_and_Space_mi_tag>()
    );
    FEDofEntity_multiIndex::index<Composite_EntType_and_Space_mi_tag>::type::iterator dit,hi_dit;
    dit = data_dofs.lower_bound(boost::make_tuple(type,space));
    hi_dit = data_dofs.upper_bound(boost::make_tuple(type,space));
    for(;dit!=hi_dit;dit++) {
      ApproximationOrder ent_order = (*dit)->get_max_order();
      int side_number = (*dit)->sideNumberPtr->side_number;
      if(side_number < 0) {
        SETERRQ(PETSC_COMM_SELF,MOFEM_DATA_INCONSISTENCY,"data inconsistency");
      }
      data[side_number].getDataOrder() = data[side_number].getDataOrder() > ent_order ? data[side_number].getDataOrder() : ent_order;
      if((*dit)->sideNumberPtr->brother_side_number!=-1) {
        if(data.size() < (unsigned int)(*dit)->sideNumberPtr->brother_side_number) {
          SETERRQ(PETSC_COMM_SELF,MOFEM_DATA_INCONSISTENCY,"data inconsistency");
        }
        data[(*dit)->sideNumberPtr->brother_side_number].getDataOrder() = data[side_number].getDataOrder();
      }
    }
  } catch (std::exception& ex) {
    std::ostringstream ss;
    ss << "thorw in method: " << ex.what() << " at line " << __LINE__ << " in file " << __FILE__;
    SETERRQ(PETSC_COMM_SELF,MOFEM_STD_EXCEPTION_THROW,ss.str().c_str());
  }
  PetscFunctionReturn(0);
}

PetscErrorCode ForcesAndSurcesCore::getEdgesDataOrder(DataForcesAndSurcesCore &data,const FieldSpace space) const {
  PetscErrorCode ierr;
  PetscFunctionBegin;
  ierr = getDataOrder(MBEDGE,space,data.dataOnEntities[MBEDGE]); CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

PetscErrorCode ForcesAndSurcesCore::getTrisDataOrder(DataForcesAndSurcesCore &data,const FieldSpace space) const {
  PetscErrorCode ierr;
  PetscFunctionBegin;
  ierr = getDataOrder(MBTRI,space,data.dataOnEntities[MBTRI]); CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

PetscErrorCode ForcesAndSurcesCore::getQuadDataOrder(DataForcesAndSurcesCore &data,const FieldSpace space) const {
  PetscErrorCode ierr;
  PetscFunctionBegin;
  ierr = getDataOrder(MBQUAD,space,data.dataOnEntities[MBQUAD]); CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

PetscErrorCode ForcesAndSurcesCore::getTetDataOrder(DataForcesAndSurcesCore &data,const FieldSpace space) const {
  PetscErrorCode ierr;
  PetscFunctionBegin;
  ierr = getDataOrder(MBTET,space,data.dataOnEntities[MBTET]); CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

PetscErrorCode ForcesAndSurcesCore::getPrismDataOrder(DataForcesAndSurcesCore &data,const FieldSpace space) const {
  PetscErrorCode ierr;
  PetscFunctionBegin;
  ierr = getDataOrder(MBPRISM,space,data.dataOnEntities[MBPRISM]); CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

PetscErrorCode ForcesAndSurcesCore::getDataOrderSpaceAndBase(
  const std::string &field_name,const EntityType type,boost::ptr_vector<DataForcesAndSurcesCore::EntData> &data
) const {
  PetscFunctionBegin;

  SideNumber_multiIndex& side_table = const_cast<SideNumber_multiIndex&>(numeredEntFiniteElementPtr->get_side_number_table());
  if(data.size() < side_table.get<2>().count(type)) {
    // prims has 9 edges, some of edges for "flat" prism are not active
    SETERRQ2(
      PETSC_COMM_SELF,
      MOFEM_DATA_INCONSISTENCY,
      "data inconsistency %d < %d",
      data.size(),side_table.get<2>().count(type)
    );
  }

  for(unsigned int side = 0;side<data.size();side++) {
    data[side].getDataOrder() = 0;
    data[side].getBase() = NOBASE;
    data[side].getSpace() = NOSPACE;
  }

  FEDofEntity_multiIndex::index<Composite_Name_And_Type_mi_tag>::type &data_dofs =
  const_cast<FEDofEntity_multiIndex::index<Composite_Name_And_Type_mi_tag>::type&>(
    numeredEntFiniteElementPtr->get_data_dofs().get<Composite_Name_And_Type_mi_tag>()
  );

  FEDofEntity_multiIndex::index<Composite_Name_And_Type_mi_tag>::type::iterator dit,hi_dit;
  dit = data_dofs.lower_bound(boost::make_tuple(field_name,type));
  hi_dit = data_dofs.upper_bound(boost::make_tuple(field_name,type));

  for(;dit!=hi_dit;dit++) {

    // std::cerr << ApproximationBaseNames[dit->get_approx_base()] << std::endl;

    int side_number = (*dit)->sideNumberPtr->side_number;
    if(data[side_number].getDataOrder()) continue;

    ApproximationOrder ent_order = (*dit)->get_max_order();
    if(side_number < 0) {
      SETERRQ(PETSC_COMM_SELF,MOFEM_DATA_INCONSISTENCY,"data inconsistency");
    }
    data[side_number].getBase() = (*dit)->get_approx_base();
    data[side_number].getSpace() = (*dit)->get_space();
    data[side_number].getDataOrder() = data[side_number].getDataOrder() > ent_order ? data[side_number].getDataOrder() : ent_order;
    if((*dit)->sideNumberPtr->brother_side_number!=-1) {
      if(data.size() < (unsigned int)(*dit)->sideNumberPtr->brother_side_number) {
        SETERRQ(PETSC_COMM_SELF,MOFEM_DATA_INCONSISTENCY,"data inconsistency");
      }
      data[(*dit)->sideNumberPtr->brother_side_number].getBase() = data[side_number].getBase();
      data[(*dit)->sideNumberPtr->brother_side_number].getSpace() = data[side_number].getSpace();
      data[(*dit)->sideNumberPtr->brother_side_number].getDataOrder() = data[side_number].getDataOrder();
    }

  }

  PetscFunctionReturn(0);
}

PetscErrorCode ForcesAndSurcesCore::getEdgesDataOrderSpaceAndBase(DataForcesAndSurcesCore &data,const std::string &field_name) const {
  PetscErrorCode ierr;
  PetscFunctionBegin;
  ierr = getDataOrderSpaceAndBase(field_name,MBEDGE,data.dataOnEntities[MBEDGE]); CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

PetscErrorCode ForcesAndSurcesCore::getTrisDataOrderSpaceAndBase(DataForcesAndSurcesCore &data,const std::string &field_name) const {
  PetscErrorCode ierr;
  PetscFunctionBegin;
  ierr = getDataOrderSpaceAndBase(field_name,MBTRI,data.dataOnEntities[MBTRI]); CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

PetscErrorCode ForcesAndSurcesCore::getQuadDataOrderSpaceAndBase(DataForcesAndSurcesCore &data,const std::string &field_name) const {
  PetscErrorCode ierr;
  PetscFunctionBegin;
  ierr = getDataOrderSpaceAndBase(field_name,MBQUAD,data.dataOnEntities[MBQUAD]); CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

PetscErrorCode ForcesAndSurcesCore::getTetDataOrderSpaceAndBase(DataForcesAndSurcesCore &data,const std::string &field_name) const {
  PetscErrorCode ierr;
  PetscFunctionBegin;
  ierr = getDataOrderSpaceAndBase(field_name,MBTET,data.dataOnEntities[MBTET]); CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

PetscErrorCode ForcesAndSurcesCore::getPrismDataOrderSpaceAndBase(DataForcesAndSurcesCore &data,const std::string &field_name) const {
  PetscErrorCode ierr;
  PetscFunctionBegin;
  ierr = getDataOrderSpaceAndBase(field_name,MBPRISM,data.dataOnEntities[MBPRISM]); CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

// ** Indices **

PetscErrorCode ForcesAndSurcesCore::getNodesIndices(
  const std::string &field_name,FENumeredDofEntity_multiIndex &dofs,VectorInt &nodes_indices,VectorInt &local_nodes_indices
) const {
  PetscErrorCode ierr;
  PetscFunctionBegin;
  FENumeredDofEntity_multiIndex::index<Composite_Name_And_Type_mi_tag>::type::iterator dit,hi_dit,it;
  dit = dofs.get<Composite_Name_And_Type_mi_tag>().lower_bound(boost::make_tuple(field_name,MBVERTEX));
  hi_dit = dofs.get<Composite_Name_And_Type_mi_tag>().upper_bound(boost::make_tuple(field_name,MBVERTEX));

  int num_nodes;
  ierr = getNumberOfNodes(num_nodes); CHKERRQ(ierr);
  int max_nb_dofs = 0;
  if(dit!=hi_dit) {
    max_nb_dofs = (*dit)->get_nb_of_coeffs()*num_nodes;
  }

  if(distance(dit,hi_dit)!=max_nb_dofs) {
    nodes_indices.resize(max_nb_dofs,false);
    local_nodes_indices.resize(max_nb_dofs,false);
    for(int dd = 0;dd<max_nb_dofs;dd++) {
      nodes_indices[dd] = -1;
      local_nodes_indices[dd] = -1;
    }
  } else {
    nodes_indices.resize(distance(dit,hi_dit),false);
    local_nodes_indices.resize(distance(dit,hi_dit),false);
  }

  for(;dit!=hi_dit;dit++) {
    int idx = (*dit)->get_petsc_gloabl_dof_idx();
    int local_idx = (*dit)->get_petsc_local_dof_idx();
    int side_number = (*dit)->sideNumberPtr->side_number;
    int pos = side_number*(*dit)->get_nb_of_coeffs()+(*dit)->get_dof_coeff_idx();
    nodes_indices[pos] = idx;
    local_nodes_indices[pos] = local_idx;
    int  brother_side_number = (*dit)->sideNumberPtr->brother_side_number;
    if(brother_side_number!=-1) {
      if(nodes_indices.size()<(unsigned int)(brother_side_number*(*dit)->get_nb_of_coeffs()+(*dit)->get_nb_of_coeffs())) {
        nodes_indices.resize(brother_side_number*(*dit)->get_nb_of_coeffs()+(*dit)->get_nb_of_coeffs());
      }
      int elem_idx = brother_side_number*(*dit)->get_nb_of_coeffs()+(*dit)->get_dof_coeff_idx();
      nodes_indices[elem_idx] = idx;
      local_nodes_indices[elem_idx] = local_idx;
    }
  }
  PetscFunctionReturn(0);
}

PetscErrorCode ForcesAndSurcesCore::getRowNodesIndices(DataForcesAndSurcesCore &data,const std::string &field_name) const {
  PetscErrorCode ierr;
  PetscFunctionBegin;
  ierr = getNodesIndices(
    field_name,const_cast<FENumeredDofEntity_multiIndex&>(numeredEntFiniteElementPtr->get_rows_dofs()),
    data.dataOnEntities[MBVERTEX][0].getIndices(),data.dataOnEntities[MBVERTEX][0].getLocalIndices()
  ); CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

PetscErrorCode ForcesAndSurcesCore::getColNodesIndices(DataForcesAndSurcesCore &data,const std::string &field_name) const {
  PetscErrorCode ierr;
  PetscFunctionBegin;
  ierr = getNodesIndices(
    field_name,const_cast<FENumeredDofEntity_multiIndex&>(numeredEntFiniteElementPtr->get_cols_dofs()),
    data.dataOnEntities[MBVERTEX][0].getIndices(),data.dataOnEntities[MBVERTEX][0].getLocalIndices()
  ); CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

PetscErrorCode ForcesAndSurcesCore::getTypeIndices(
  const std::string &field_name,FENumeredDofEntity_multiIndex &dofs,EntityType type,int side_number,VectorInt &indices,VectorInt &local_indices
) const {
  PetscErrorCode ierr;
  PetscFunctionBegin;
  FENumeredDofEntity_multiIndex::index<Composite_Name_Type_And_Side_Number_mi_tag>::type::iterator dit,hi_dit;
  dit = dofs.get<Composite_Name_Type_And_Side_Number_mi_tag>().lower_bound(boost::make_tuple(field_name,type,side_number));
  hi_dit = dofs.get<Composite_Name_Type_And_Side_Number_mi_tag>().upper_bound(boost::make_tuple(field_name,type,side_number));
  if(dit!=hi_dit) {
    indices.resize((*dit)->get_nb_dofs_on_ent(),false);
    local_indices.resize((*dit)->get_nb_dofs_on_ent(),false);
    for(;dit!=hi_dit;dit++) {
      int idx = (*dit)->get_petsc_gloabl_dof_idx();
      int elemem_idx = (*dit)->get_EntDofIdx();
      indices[elemem_idx] = idx;
      int local_idx = (*dit)->get_petsc_local_dof_idx();
      local_indices[elemem_idx] = local_idx;
    }
  } else {
    indices.resize(0);
    local_indices.resize(0);
  }
  PetscFunctionReturn(0);
}

PetscErrorCode ForcesAndSurcesCore::getTypeIndices(
  const std::string &field_name,FENumeredDofEntity_multiIndex &dofs,EntityType type,
  boost::ptr_vector<DataForcesAndSurcesCore::EntData> &data
) const {
  PetscErrorCode ierr;
  PetscFunctionBegin;
  SideNumber_multiIndex& side_table = const_cast<SideNumber_multiIndex&>(numeredEntFiniteElementPtr->get_side_number_table());
  SideNumber_multiIndex::nth_index<2>::type::iterator siit = side_table.get<2>().lower_bound(type);
  SideNumber_multiIndex::nth_index<2>::type::iterator hi_siit = side_table.get<2>().upper_bound(type);
  for(;siit!=hi_siit;siit++) {
    ierr = getTypeIndices(
      field_name,dofs,type,siit->get()->side_number,data[siit->get()->side_number].getIndices(),data[siit->get()->side_number].getLocalIndices()
    ); CHKERRQ(ierr);
    if(siit->get()->brother_side_number!=-1) {
      ierr = getTypeIndices(
        field_name,dofs,type,siit->get()->side_number,data[siit->get()->brother_side_number].getIndices(),data[siit->get()->brother_side_number].getLocalIndices()
      ); CHKERRQ(ierr);
    }
  }
  PetscFunctionReturn(0);
}

PetscErrorCode ForcesAndSurcesCore::getEdgesRowIndices(DataForcesAndSurcesCore &data,const std::string &field_name) const {
  PetscErrorCode ierr;
  PetscFunctionBegin;
  ierr = getTypeIndices(
    field_name,const_cast<FENumeredDofEntity_multiIndex&>(numeredEntFiniteElementPtr->get_rows_dofs()),MBEDGE,data.dataOnEntities[MBEDGE]
  ); CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

PetscErrorCode ForcesAndSurcesCore::getEdgesColIndices(DataForcesAndSurcesCore &data,const std::string &field_name) const {
  PetscErrorCode ierr;
  PetscFunctionBegin;
  ierr = getTypeIndices(
    field_name,const_cast<FENumeredDofEntity_multiIndex&>(numeredEntFiniteElementPtr->get_cols_dofs()),MBEDGE,data.dataOnEntities[MBEDGE]
  ); CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

PetscErrorCode ForcesAndSurcesCore::getTrisRowIndices(
  DataForcesAndSurcesCore &data,const std::string &field_name
) const {
  PetscErrorCode ierr;
  PetscFunctionBegin;
  ierr = getTypeIndices(
    field_name,const_cast<FENumeredDofEntity_multiIndex&>(numeredEntFiniteElementPtr->get_rows_dofs()),MBTRI,data.dataOnEntities[MBTRI]
  ); CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

PetscErrorCode ForcesAndSurcesCore::getTrisColIndices(DataForcesAndSurcesCore &data,const std::string &field_name) const {
  PetscErrorCode ierr;
  PetscFunctionBegin;
  ierr = getTypeIndices(
    field_name,const_cast<FENumeredDofEntity_multiIndex&>(numeredEntFiniteElementPtr->get_cols_dofs()),MBTRI,data.dataOnEntities[MBTRI]
  ); CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

PetscErrorCode ForcesAndSurcesCore::getTetsRowIndices(DataForcesAndSurcesCore &data,const std::string &field_name) const {
  PetscErrorCode ierr;
  PetscFunctionBegin;
  if(data.dataOnEntities[MBTET].size() == 0) {
    SETERRQ(PETSC_COMM_SELF,MOFEM_DATA_INCONSISTENCY,"data inconsistency");
  }
  ierr = getTypeIndices(
    field_name,const_cast<FENumeredDofEntity_multiIndex&>(numeredEntFiniteElementPtr->get_rows_dofs()),MBTET,0,
    data.dataOnEntities[MBTET][0].getIndices(),data.dataOnEntities[MBTET][0].getLocalIndices()
  ); CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

PetscErrorCode ForcesAndSurcesCore::getTetsColIndices(DataForcesAndSurcesCore &data,const std::string &field_name) const {
  PetscErrorCode ierr;
  PetscFunctionBegin;
  if(data.dataOnEntities[MBTET].size() == 0) {
    SETERRQ(PETSC_COMM_SELF,MOFEM_DATA_INCONSISTENCY,"data inconsistency");
  }
  ierr = getTypeIndices(
    field_name,const_cast<FENumeredDofEntity_multiIndex&>(
      numeredEntFiniteElementPtr->get_cols_dofs()),MBTET,0,data.dataOnEntities[MBTET][0].getIndices(),data.dataOnEntities[MBTET][0].getLocalIndices()
  ); CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

PetscErrorCode ForcesAndSurcesCore::getQuadRowIndices(
  DataForcesAndSurcesCore &data,const std::string &field_name
) const {
  PetscErrorCode ierr;
  PetscFunctionBegin;
  ierr = getTypeIndices(
    field_name,const_cast<FENumeredDofEntity_multiIndex&>(numeredEntFiniteElementPtr->get_rows_dofs()),MBQUAD,data.dataOnEntities[MBQUAD]
  ); CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

PetscErrorCode ForcesAndSurcesCore::getQuadColIndices(DataForcesAndSurcesCore &data,const std::string &field_name) const {
  PetscErrorCode ierr;
  PetscFunctionBegin;
  ierr = getTypeIndices(
    field_name,const_cast<FENumeredDofEntity_multiIndex&>(numeredEntFiniteElementPtr->get_cols_dofs()),MBQUAD,data.dataOnEntities[MBQUAD]
  ); CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

PetscErrorCode ForcesAndSurcesCore::getPrismRowIndices(
  DataForcesAndSurcesCore &data,const std::string &field_name
) const {
  PetscErrorCode ierr;
  PetscFunctionBegin;
  ierr = getTypeIndices(
    field_name,const_cast<FENumeredDofEntity_multiIndex&>(numeredEntFiniteElementPtr->get_rows_dofs()),MBPRISM,data.dataOnEntities[MBPRISM]
  ); CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

PetscErrorCode ForcesAndSurcesCore::getPrismColIndices(DataForcesAndSurcesCore &data,const std::string &field_name) const {
  PetscErrorCode ierr;
  PetscFunctionBegin;
  ierr = getTypeIndices(
    field_name,const_cast<FENumeredDofEntity_multiIndex&>(numeredEntFiniteElementPtr->get_cols_dofs()),MBPRISM,data.dataOnEntities[MBPRISM]
  ); CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

PetscErrorCode ForcesAndSurcesCore::getNoFieldIndices(
  const std::string &field_name,FENumeredDofEntity_multiIndex &dofs,VectorInt &indices
) const {
  PetscErrorCode ierr;
  PetscFunctionBegin;
  FENumeredDofEntity_multiIndex::index<FieldName_mi_tag>::type::iterator dit,hi_dit;
  dit = dofs.get<FieldName_mi_tag>().lower_bound(field_name);
  hi_dit = dofs.get<FieldName_mi_tag>().upper_bound(field_name);
  indices.resize(distance(dit,hi_dit));
  for(;dit!=hi_dit;dit++) {
    int idx = (*dit)->get_petsc_gloabl_dof_idx();
    indices[(*dit)->get_dof_coeff_idx()] = idx;
  }
  PetscFunctionReturn(0);
}

PetscErrorCode ForcesAndSurcesCore::getNoFieldRowIndices(DataForcesAndSurcesCore &data,const std::string &field_name) const {
  PetscErrorCode ierr;
  PetscFunctionBegin;
  //EntityType fe_type = numeredEntFiniteElementPtr->get_ent_type();
  if(data.dataOnEntities[MBENTITYSET].size() == 0) {
    SETERRQ(PETSC_COMM_SELF,MOFEM_DATA_INCONSISTENCY,"data inconsistency");
  }
  ierr = getNoFieldIndices(
    field_name,const_cast<FENumeredDofEntity_multiIndex&>(numeredEntFiniteElementPtr->get_rows_dofs()),data.dataOnEntities[MBENTITYSET][0].getIndices()
  ); CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

PetscErrorCode ForcesAndSurcesCore::getNoFieldColIndices(DataForcesAndSurcesCore &data,const std::string &field_name) const {
  PetscErrorCode ierr;
  PetscFunctionBegin;
  if(data.dataOnEntities[MBENTITYSET].size() == 0) {
    SETERRQ(PETSC_COMM_SELF,MOFEM_DATA_INCONSISTENCY,"data inconsistency");
  }
  ierr = getNoFieldIndices(
    field_name,const_cast<FENumeredDofEntity_multiIndex&>(numeredEntFiniteElementPtr->get_cols_dofs()),data.dataOnEntities[MBENTITYSET][0].getIndices()
  ); CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

// ** Indices from problem **

PetscErrorCode ForcesAndSurcesCore::getProblemNodesIndices(
  const std::string &field_name,
  const NumeredDofEntity_multiIndex &dofs,
  VectorInt &nodes_indices
) const {
  PetscFunctionBegin;

  nodes_indices.resize(0);

  SideNumber_multiIndex& side_table = const_cast<SideNumber_multiIndex&>(numeredEntFiniteElementPtr->get_side_number_table());
  SideNumber_multiIndex::nth_index<1>::type::iterator siit = side_table.get<1>().lower_bound(boost::make_tuple(MBVERTEX,0));
  SideNumber_multiIndex::nth_index<1>::type::iterator hi_siit = side_table.get<1>().lower_bound(boost::make_tuple(MBVERTEX,10000));

  int nn = 0;
  for(;siit!=hi_siit;siit++,nn++) {

    if(siit->get()->side_number == -1) continue;

    const EntityHandle ent = siit->get()->ent;
    NumeredDofEntity_multiIndex::index<Composite_Name_And_Ent_And_EndDofIdx_mi_tag>::type::iterator dit,hi_dit;
    dit = dofs.get<Composite_Name_And_Ent_And_EndDofIdx_mi_tag>().lower_bound(boost::make_tuple(field_name,ent,0));
    hi_dit = dofs.get<Composite_Name_And_Ent_And_EndDofIdx_mi_tag>().upper_bound(boost::make_tuple(field_name,ent,10000));  /// very large number

    if(dit!=hi_dit) {

      if(!nn) {
        nodes_indices.resize((*dit)->get_nb_of_coeffs()*distance(siit,hi_siit));
      }
      for(;dit!=hi_dit;dit++) {
        nodes_indices[siit->get()->side_number*(*dit)->get_nb_of_coeffs()+(*dit)->get_dof_coeff_idx()] = (*dit)->get_petsc_gloabl_dof_idx();
      }

    }

  }

  PetscFunctionReturn(0);

}

PetscErrorCode ForcesAndSurcesCore::getProblemTypeIndices(
  const std::string &field_name,const NumeredDofEntity_multiIndex &dofs,
  EntityType type,int side_number,VectorInt &indices
) const {
  PetscFunctionBegin;

  indices.resize(0);

  SideNumber_multiIndex& side_table = const_cast<SideNumber_multiIndex&>(numeredEntFiniteElementPtr->get_side_number_table());
  SideNumber_multiIndex::nth_index<1>::type::iterator siit = side_table.get<1>().lower_bound(boost::make_tuple(type,side_number));
  SideNumber_multiIndex::nth_index<1>::type::iterator hi_siit = side_table.get<1>().upper_bound(boost::make_tuple(type,side_number));

  for(;siit!=hi_siit;siit++) {

    if(siit->get()->side_number == -1) continue;

    const EntityHandle ent = siit->get()->ent;
    NumeredDofEntity_multiIndex::index<Composite_Name_And_Ent_And_EndDofIdx_mi_tag>::type::iterator dit,hi_dit;
    dit = dofs.get<Composite_Name_And_Ent_And_EndDofIdx_mi_tag>().lower_bound(boost::make_tuple(field_name,ent,0));
    hi_dit = dofs.get<Composite_Name_And_Ent_And_EndDofIdx_mi_tag>().upper_bound(boost::make_tuple(field_name,ent,10000));  /// very large number

    indices.resize(distance(dit,hi_dit));
    for(;dit!=hi_dit;dit++) {

      indices[(*dit)->get_EntDofIdx()] = (*dit)->get_petsc_gloabl_dof_idx();

    }

  }

  PetscFunctionReturn(0);
}

PetscErrorCode ForcesAndSurcesCore::getProblemNodesRowIndices(const std::string &field_name,VectorInt &nodes_indices) const {
  PetscFunctionBegin;
  PetscErrorCode ierr;
  ierr = getProblemNodesIndices(field_name,*(problemPtr->numered_dofs_rows),nodes_indices); CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

PetscErrorCode ForcesAndSurcesCore::getProblemTypeRowIndices(const std::string &field_name,EntityType type,int side_number,VectorInt &indices) const {
  PetscFunctionBegin;
  PetscErrorCode ierr;
  ierr = getProblemTypeIndices(field_name,*(problemPtr->numered_dofs_rows),type,side_number,indices); CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

PetscErrorCode ForcesAndSurcesCore::getProblemNodesColIndices(const std::string &field_name,VectorInt &nodes_indices) const {
  PetscFunctionBegin;
  PetscErrorCode ierr;
  ierr = getProblemNodesIndices(field_name,*(problemPtr->numered_dofs_cols),nodes_indices); CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

PetscErrorCode ForcesAndSurcesCore::getProblemTypeColIndices(const std::string &field_name,EntityType type,int side_number,VectorInt &indices) const {
  PetscFunctionBegin;
  PetscErrorCode ierr;
  ierr = getProblemTypeIndices(field_name,*(problemPtr->numered_dofs_cols),type,side_number,indices); CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

// ** Data **

PetscErrorCode ForcesAndSurcesCore::getNodesFieldData(
  const std::string &field_name,
  FEDofEntity_multiIndex &dofs,
  VectorDouble &nodes_data,
  VectorDofs &nodes_dofs,
  FieldSpace &space,
  FieldApproximationBase &base
) const {
  PetscErrorCode ierr;
  PetscFunctionBegin;
  try {
    FEDofEntity_multiIndex::index<Composite_Name_And_Type_mi_tag>::type::iterator dit,hi_dit,it;
    dit = dofs.get<Composite_Name_And_Type_mi_tag>().lower_bound(boost::make_tuple(field_name,MBVERTEX));
    hi_dit = dofs.get<Composite_Name_And_Type_mi_tag>().upper_bound(boost::make_tuple(field_name,MBVERTEX));

    int num_nodes;
    ierr = getNumberOfNodes(num_nodes); CHKERRQ(ierr);
    int max_nb_dofs = 0;
    if(dit!=hi_dit) {
      max_nb_dofs = (*dit)->get_nb_of_coeffs()*num_nodes;
    }

    if(distance(dit,hi_dit)!=max_nb_dofs) {
      nodes_data.resize(max_nb_dofs,false);
      nodes_data.clear();
      nodes_dofs.resize(max_nb_dofs,false);
    } else {
      int size = distance(dit,hi_dit);
      nodes_data.resize(size,false);
      nodes_dofs.resize(size,false);
    }

    if(dit!=hi_dit) {
      space = (*dit)->get_space();
      base = (*dit)->get_approx_base();
    }

    for(;dit!=hi_dit;dit++) {
      FieldData val = (*dit)->get_FieldData();
      int side_number = (*dit)->sideNumberPtr->side_number;
      if(side_number == -1) {
        SETERRQ(PETSC_COMM_SELF,MOFEM_DATA_INCONSISTENCY,"data inconsistency");
      }
      int pos = side_number*(*dit)->get_nb_of_coeffs()+(*dit)->get_dof_coeff_idx();
      nodes_data[pos] = val;
      nodes_dofs[pos] = &*(*dit);
      int  brother_side_number = (*dit)->sideNumberPtr->brother_side_number;
      if(brother_side_number!=-1) {
        if(nodes_data.size()<(unsigned int)(brother_side_number*(*dit)->get_nb_of_coeffs()+(*dit)->get_nb_of_coeffs())) {
          nodes_data.resize(brother_side_number*(*dit)->get_nb_of_coeffs()+(*dit)->get_nb_of_coeffs());
          nodes_dofs.resize(brother_side_number*(*dit)->get_nb_of_coeffs()+(*dit)->get_nb_of_coeffs());
        }
        int brother_pos = brother_side_number*(*dit)->get_nb_of_coeffs()+(*dit)->get_dof_coeff_idx();
        nodes_data[brother_pos] = val;
        nodes_dofs[brother_pos] = &*(*dit);
      }
    }

  } catch (std::exception& ex) {
    std::ostringstream ss;
    ss << "thorw in method: " << ex.what() << " at line " << __LINE__ << " in file " << __FILE__;
    SETERRQ(PETSC_COMM_SELF,MOFEM_STD_EXCEPTION_THROW,ss.str().c_str());
  }
  PetscFunctionReturn(0);
}

PetscErrorCode ForcesAndSurcesCore::getNodesFieldData(DataForcesAndSurcesCore &data,const std::string &field_name) const {
  PetscFunctionBegin;
  PetscErrorCode ierr;
  ierr = getNodesFieldData(
    field_name,
    const_cast<FEDofEntity_multiIndex&>(numeredEntFiniteElementPtr->get_data_dofs()),
    data.dataOnEntities[MBVERTEX][0].getFieldData(),
    data.dataOnEntities[MBVERTEX][0].getFieldDofs(),
    data.dataOnEntities[MBVERTEX][0].getSpace(),
    data.dataOnEntities[MBVERTEX][0].getBase()
  ); CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

PetscErrorCode ForcesAndSurcesCore::getTypeFieldData(
  const std::string &field_name,
  FEDofEntity_multiIndex &dofs,
  EntityType type,
  int side_number,
  VectorDouble &ent_field_data,
  VectorDofs &ent_field_dofs
) const {
  PetscErrorCode ierr;
  PetscFunctionBegin;
  FEDofEntity_multiIndex::index<Composite_Name_Type_And_Side_Number_mi_tag>::type::iterator dit,hi_dit;
  dit = dofs.get<Composite_Name_Type_And_Side_Number_mi_tag>().lower_bound(boost::make_tuple(field_name,type,side_number));
  hi_dit = dofs.get<Composite_Name_Type_And_Side_Number_mi_tag>().upper_bound(boost::make_tuple(field_name,type,side_number));
  if(dit!=hi_dit) {
    ent_field_data.resize((*dit)->get_nb_dofs_on_ent(),false);
    ent_field_dofs.resize((*dit)->get_nb_dofs_on_ent(),false);
    for(;dit!=hi_dit;dit++) {
      FieldData val = (*dit)->get_FieldData();
      int idx = (*dit)->get_EntDofIdx();
      ent_field_data[idx] = val;
      ent_field_dofs[idx] = &*(*dit);
    }
  } else {
    ent_field_data.resize(0);
    ent_field_dofs.resize(0);
  }
  PetscFunctionReturn(0);
}

PetscErrorCode ForcesAndSurcesCore::getTypeFieldData(
  const std::string &field_name,
  FEDofEntity_multiIndex &dofs,
  EntityType type,
  boost::ptr_vector<DataForcesAndSurcesCore::EntData> &data
) const {
  PetscErrorCode ierr;
  PetscFunctionBegin;
  SideNumber_multiIndex& side_table = const_cast<SideNumber_multiIndex&>(numeredEntFiniteElementPtr->get_side_number_table());
  if(data.size() < side_table.get<2>().count(type)) {
    SETERRQ(PETSC_COMM_SELF,MOFEM_DATA_INCONSISTENCY,"data inconsistency");
  }
  SideNumber_multiIndex::nth_index<2>::type::iterator siit = side_table.get<2>().lower_bound(type);
  SideNumber_multiIndex::nth_index<2>::type::iterator hi_siit = side_table.get<2>().upper_bound(type);
  for(;siit!=hi_siit;siit++) {
    ierr = getTypeFieldData(
      field_name,dofs,type,
      siit->get()->side_number,
      data[siit->get()->side_number].getFieldData(),
      data[siit->get()->side_number].getFieldDofs()
    ); CHKERRQ(ierr);
    if(siit->get()->brother_side_number!=-1) {
      if(data.size() < (unsigned int)siit->get()->brother_side_number) {
        SETERRQ(PETSC_COMM_SELF,MOFEM_DATA_INCONSISTENCY,"data inconsistency");
      }
      ierr = getTypeFieldData(
        field_name,dofs,type,siit->get()->side_number,
        data[siit->get()->brother_side_number].getFieldData(),
        data[siit->get()->brother_side_number].getFieldDofs()
      ); CHKERRQ(ierr);
    }
  }
  PetscFunctionReturn(0);
}

PetscErrorCode ForcesAndSurcesCore::getNoFieldFieldData(
  const std::string &field_name,
  FEDofEntity_multiIndex &dofs,
  VectorDouble &ent_field_data,
  VectorDofs &ent_field_dofs
) const {
  PetscErrorCode ierr;
  PetscFunctionBegin;
  FEDofEntity_multiIndex::index<FieldName_mi_tag>::type::iterator dit,hi_dit;
  dit = dofs.get<FieldName_mi_tag>().lower_bound(field_name);
  hi_dit = dofs.get<FieldName_mi_tag>().upper_bound(field_name);
  int size = distance(dit,hi_dit);
  ent_field_data.resize(size,false);
  ent_field_dofs.resize(size,false);
  for(;dit!=hi_dit;dit++) {
    int idx = (*dit)->get_dof_coeff_idx();
    ent_field_data[idx] = (*dit)->get_FieldData();
    ent_field_dofs[idx] = &*(*dit);
  }
  PetscFunctionReturn(0);
}

PetscErrorCode ForcesAndSurcesCore::getNoFieldFieldData(
  DataForcesAndSurcesCore &data,const std::string &field_name
) const {
  PetscErrorCode ierr;
  PetscFunctionBegin;
  if(data.dataOnEntities[MBENTITYSET].size() == 0) {
    SETERRQ(PETSC_COMM_SELF,MOFEM_DATA_INCONSISTENCY,"data inconsistency");
  }
  ierr = getNoFieldFieldData(
    field_name,const_cast<FEDofEntity_multiIndex&>(numeredEntFiniteElementPtr->get_data_dofs()),
    data.dataOnEntities[MBENTITYSET][0].getFieldData(),
    data.dataOnEntities[MBENTITYSET][0].getFieldDofs()
  ); CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

PetscErrorCode ForcesAndSurcesCore::getEdgesFieldData(DataForcesAndSurcesCore &data,const std::string &field_name) const {
  PetscErrorCode ierr;
  PetscFunctionBegin;
  ierr = getTypeFieldData(
    field_name,const_cast<FEDofEntity_multiIndex&>(numeredEntFiniteElementPtr->get_data_dofs()),MBEDGE,data.dataOnEntities[MBEDGE]
  ); CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

PetscErrorCode ForcesAndSurcesCore::getTrisFieldData(
  DataForcesAndSurcesCore &data,const std::string &field_name
) const {
  PetscErrorCode ierr;
  PetscFunctionBegin;
  ierr = getTypeFieldData(
    field_name,const_cast<FEDofEntity_multiIndex&>(numeredEntFiniteElementPtr->get_data_dofs()),MBTRI,data.dataOnEntities[MBTRI]
  ); CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

PetscErrorCode ForcesAndSurcesCore::getQuadFieldData(
  DataForcesAndSurcesCore &data,const std::string &field_name
) const {
  PetscErrorCode ierr;
  PetscFunctionBegin;
  ierr = getTypeFieldData(
    field_name,const_cast<FEDofEntity_multiIndex&>(numeredEntFiniteElementPtr->get_data_dofs()),MBQUAD,data.dataOnEntities[MBQUAD]
  ); CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

PetscErrorCode ForcesAndSurcesCore::getTetsFieldData(DataForcesAndSurcesCore &data,const std::string &field_name) const {
  PetscErrorCode ierr;
  PetscFunctionBegin;
  if(data.dataOnEntities[MBTET].size() == 0) {
    SETERRQ(PETSC_COMM_SELF,MOFEM_DATA_INCONSISTENCY,"data inconsistency");
  }
  ierr = getTypeFieldData(
    field_name,const_cast<FEDofEntity_multiIndex&>(numeredEntFiniteElementPtr->get_data_dofs()),
    MBTET,
    0,
    data.dataOnEntities[MBTET][0].getFieldData(),
    data.dataOnEntities[MBTET][0].getFieldDofs()
  ); CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

PetscErrorCode ForcesAndSurcesCore::getPrismFieldData(DataForcesAndSurcesCore &data,const std::string &field_name) const {
  PetscErrorCode ierr;
  PetscFunctionBegin;
  if(data.dataOnEntities[MBPRISM].size() == 0) {
    SETERRQ(PETSC_COMM_SELF,MOFEM_DATA_INCONSISTENCY,"data inconsistency");
  }
  ierr = getTypeFieldData(
    field_name,
    const_cast<FEDofEntity_multiIndex&>(numeredEntFiniteElementPtr->get_data_dofs()),
    MBPRISM,
    0,
    data.dataOnEntities[MBPRISM][0].getFieldData(),
    data.dataOnEntities[MBPRISM][0].getFieldDofs()
  ); CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

// ** Face **

PetscErrorCode ForcesAndSurcesCore::getFaceTriNodes(DataForcesAndSurcesCore &data) const {
  PetscErrorCode ierr;
  PetscFunctionBegin;
  //PetscAttachDebugger();
  data.facesNodes.resize(4,3,false);
  SideNumber_multiIndex& side_table = const_cast<SideNumber_multiIndex&>(numeredEntFiniteElementPtr->get_side_number_table());
  SideNumber_multiIndex::nth_index<1>::type::iterator siit = side_table.get<1>().lower_bound(boost::make_tuple(MBTRI,0));
  SideNumber_multiIndex::nth_index<1>::type::iterator hi_siit = side_table.get<1>().upper_bound(boost::make_tuple(MBTRI,4));
  const int cannonical_face_sense_p1[4][3] = { {0,1,3}, {1,2,3}, {0,3,2}/**/, {0,2,1}/**/ }; //secon index is offset (positive sense)
  const int cannonical_face_sense_m1[4][3] = { {0,3,1}, {1,3,2}, {0,2,3}, {0,1,2} }; //second index is offset (negative sense
  for(;siit!=hi_siit;siit++) {
    const boost::shared_ptr<SideNumber> side = *siit;
    int face_conn[3] = {-1,-1,-1};
    if(side->offset == 0) {
      face_conn[0] = side->sense == 1 ? cannonical_face_sense_p1[side->side_number][0] : cannonical_face_sense_m1[side->side_number][0];
      face_conn[1] = side->sense == 1 ? cannonical_face_sense_p1[side->side_number][1] : cannonical_face_sense_m1[side->side_number][1];
      face_conn[2] = side->sense == 1 ? cannonical_face_sense_p1[side->side_number][2] : cannonical_face_sense_m1[side->side_number][2];
    }
    if(side->offset == 1) {
      face_conn[0] = side->sense == 1 ? cannonical_face_sense_p1[side->side_number][1] : cannonical_face_sense_m1[side->side_number][2]/**/;
      face_conn[1] = side->sense == 1 ? cannonical_face_sense_p1[side->side_number][2] : cannonical_face_sense_m1[side->side_number][0];
      face_conn[2] = side->sense == 1 ? cannonical_face_sense_p1[side->side_number][0] : cannonical_face_sense_m1[side->side_number][1];
    }
    if(side->offset == 2) {
      face_conn[0] = side->sense == 1 ? cannonical_face_sense_p1[side->side_number][2] : cannonical_face_sense_m1[side->side_number][1]/**/;
      face_conn[1] = side->sense == 1 ? cannonical_face_sense_p1[side->side_number][0] : cannonical_face_sense_m1[side->side_number][2];
      face_conn[2] = side->sense == 1 ? cannonical_face_sense_p1[side->side_number][1] : cannonical_face_sense_m1[side->side_number][0];
    }
    for(int nn = 0;nn<3;nn++) data.facesNodes(side->side_number,nn) = face_conn[nn];
    {
      const EntityHandle *conn_tet;
      int num_nodes_tet;
      EntityHandle ent = numeredEntFiniteElementPtr->get_ent();
      rval = mField.get_moab().get_connectivity(ent,conn_tet,num_nodes_tet,true); CHKERRQ_MOAB(rval);
      if(num_nodes_tet != 4) SETERRQ(PETSC_COMM_SELF,MOFEM_DATA_INCONSISTENCY,"data inconsistency");
      int num_nodes_face;
      const EntityHandle *conn_face;
      rval = mField.get_moab().get_connectivity(side->ent,conn_face,num_nodes_face,true); CHKERRQ_MOAB(rval);
      if(num_nodes_face != 3) SETERRQ(PETSC_COMM_SELF,1,"data inconsistency");
      if(conn_face[0] != conn_tet[data.facesNodes(side->side_number,0)])
        SETERRQ(PETSC_COMM_SELF,MOFEM_DATA_INCONSISTENCY,"data inconsistency");
      if(conn_face[1] != conn_tet[data.facesNodes(side->side_number,1)])
        SETERRQ(PETSC_COMM_SELF,MOFEM_DATA_INCONSISTENCY,"data inconsistency");
      if(conn_face[2] != conn_tet[data.facesNodes(side->side_number,2)])
        SETERRQ(PETSC_COMM_SELF,MOFEM_DATA_INCONSISTENCY,"data inconsistency");
    }
  }
  PetscFunctionReturn(0);
}

// ** Space and Base **

PetscErrorCode ForcesAndSurcesCore::getSpacesAndBaseOnEntities(DataForcesAndSurcesCore &data) const {
  PetscErrorCode ierr;
  PetscFunctionBegin;
  try {
    for(_IT_GET_FEDATA_DOFS_FOR_LOOP_(this,dof)) {
      // std::cerr << *dof << std::endl;
      // std::cerr << dof->get_space() << " " << data.sPace.size() << std::endl;
      // std::cerr << dof->get_approx_base() << " " << data.bAse.size() << std::endl;
      data.sPace.set((*dof)->get_space());
      data.bAse.set((*dof)->get_approx_base());
      data.spacesOnEntities[(*dof)->get_ent_type()].set((*dof)->get_space());
      data.basesOnEntities[(*dof)->get_ent_type()].set((*dof)->get_approx_base());
      // std::cerr << "approx base " << ApproximationBaseNames[dof->get_approx_base()] << " " << data.basesOnEntities[dof->get_ent_type()] << std::endl;
    }
  } catch (std::exception& ex) {
    std::ostringstream ss;
    ss << "thorw in method: " << ex.what() << " at line " << __LINE__ << " in file " << __FILE__;
    SETERRQ(PETSC_COMM_SELF,MOFEM_STD_EXCEPTION_THROW,ss.str().c_str());
  }
  PetscFunctionReturn(0);
}

// **** Data Operator ****

static PetscErrorCode get_porblem_row_indices(
  const ForcesAndSurcesCore *fe_ptr,const EntityType type,const int side,const std::string field_name,VectorInt& indices) {
  PetscFunctionBegin;

  PetscErrorCode ierr;

  switch(type) {
    case MBVERTEX:
      ierr = fe_ptr->getProblemNodesRowIndices(field_name,indices); CHKERRQ(ierr);
    break;
    default:
      ierr = fe_ptr->getProblemTypeRowIndices(field_name,type,side,indices); CHKERRQ(ierr);
  }

  PetscFunctionReturn(0);
}

static PetscErrorCode get_porblem_col_indices(
  const ForcesAndSurcesCore *fe_ptr,const EntityType type,const int side,const std::string field_name,VectorInt& indices) {
  PetscFunctionBegin;;

  PetscErrorCode ierr;

  switch(type) {
    case MBVERTEX:
      ierr = fe_ptr->getProblemNodesColIndices(field_name,indices); CHKERRQ(ierr);
    break;
    default:
      ierr = fe_ptr->getProblemTypeColIndices(field_name,type,side,indices); CHKERRQ(ierr);
  }

  PetscFunctionReturn(0);
}

PetscErrorCode ForcesAndSurcesCore::UserDataOperator::getPorblemRowIndices(
  const std::string field_name,const EntityType type,const int side,VectorInt& indices
) const {
  PetscFunctionBegin;

  PetscErrorCode ierr;

  if(ptrFE == NULL) {
    SETERRQ(PETSC_COMM_SELF,MOFEM_DATA_INCONSISTENCY,"data inconsistency");
  }
  ierr = get_porblem_row_indices(ptrFE,type,side,field_name,indices); CHKERRQ(ierr);

  PetscFunctionReturn(0);
}

PetscErrorCode ForcesAndSurcesCore::UserDataOperator::getPorblemColIndices(
  const std::string field_name,const EntityType type,const int side,VectorInt& indices
) const {
  PetscFunctionBegin;

  PetscErrorCode ierr;

  if(ptrFE == NULL) {
    SETERRQ(PETSC_COMM_SELF,MOFEM_DATA_INCONSISTENCY,"data inconsistency");
  }
  ierr = get_porblem_col_indices(ptrFE,type,side,field_name,indices); CHKERRQ(ierr);

  PetscFunctionReturn(0);
}

}
