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

PetscErrorCode ForcesAndSurcesCore::getNumberOfNodes(int &num_nodes) {
  PetscFunctionBegin;

  EntityHandle ent = fePtr->get_ent();
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

PetscErrorCode ForcesAndSurcesCore::getSense(EntityType type,boost::ptr_vector<DataForcesAndSurcesCore::EntData> &data) {
  PetscFunctionBegin;
  try {
    SideNumber_multiIndex& side_table = const_cast<SideNumber_multiIndex&>(fePtr->get_side_number_table());
    if(data.size() < side_table.get<2>().count(type)) {
      // prims has 9 edges, some of edges for "flat" prism are not active
      SETERRQ2(PETSC_COMM_SELF,MOFEM_DATA_INCONSISTENCY,"data inconsistency %u < %u",data.size(),side_table.get<2>().count(type));
    }
    SideNumber_multiIndex::nth_index<2>::type::iterator siit = side_table.get<2>().lower_bound(type);
    SideNumber_multiIndex::nth_index<2>::type::iterator hi_siit = side_table.get<2>().upper_bound(type);
    for(;siit!=hi_siit;siit++) {
      data[siit->side_number].getSense() = siit->sense;
      if(siit->brother_side_number!=-1) {
        if(data.size() < (unsigned)siit->brother_side_number) {
          SETERRQ(PETSC_COMM_SELF,MOFEM_DATA_INCONSISTENCY,"data inconsistency");
        }
        data[siit->brother_side_number].getSense() = siit->sense;
      }
    }
  } catch (exception& ex) {
    ostringstream ss;
    ss << "thorw in method: " << ex.what() << " at line " << __LINE__ << " in file " << __FILE__;
    SETERRQ(PETSC_COMM_SELF,MOFEM_STD_EXCEPTION_THROW,ss.str().c_str());
  }

  PetscFunctionReturn(0);
}

PetscErrorCode ForcesAndSurcesCore::getEdgesSense(DataForcesAndSurcesCore &data) {
  PetscFunctionBegin;
  ierr = getSense(MBEDGE,data.dataOnEntities[MBEDGE]); CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

PetscErrorCode ForcesAndSurcesCore::getTrisSense(DataForcesAndSurcesCore &data) {
  PetscFunctionBegin;
  ierr = getSense(MBTRI,data.dataOnEntities[MBTRI]); CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

PetscErrorCode ForcesAndSurcesCore::getQuadSense(DataForcesAndSurcesCore &data) {
  PetscFunctionBegin;
  ierr = getSense(MBQUAD,data.dataOnEntities[MBQUAD]); CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

// ** Order **

template<typename DOFMULTIINDEX>
static int get_max_order(
  DOFMULTIINDEX &dof_multi_index
) {
  int max_order = 0;
  typename DOFMULTIINDEX::iterator dit,hi_dit;
  dit = dof_multi_index.begin();
  hi_dit = dof_multi_index.end();
  for(;dit!=hi_dit;dit++) {
    int dit_max_order = dit->get_max_order();
    max_order = (max_order>dit_max_order) ? max_order : dit_max_order;
  }
  return max_order;
}

int ForcesAndSurcesCore::getMaxDataOrder() {
  return get_max_order(fePtr->get_data_dofs());
}

int ForcesAndSurcesCore::getMaxRowOrder() {
  return get_max_order(fePtr->get_rows_dofs());
}

int ForcesAndSurcesCore::getMaxColOrder() {
  return get_max_order(fePtr->get_cols_dofs());
}

PetscErrorCode ForcesAndSurcesCore::getDataOrder(const EntityType type,const FieldSpace space,boost::ptr_vector<DataForcesAndSurcesCore::EntData> &data) {
  PetscFunctionBegin;
  try {
    SideNumber_multiIndex& side_table = const_cast<SideNumber_multiIndex&>(fePtr->get_side_number_table());
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
    FEDofMoFEMEntity_multiIndex::index<Composite_EntType_and_Space_mi_tag>::type &data_dofs =
    const_cast<FEDofMoFEMEntity_multiIndex::index<Composite_EntType_and_Space_mi_tag>::type&>(
      fePtr->get_data_dofs().get<Composite_EntType_and_Space_mi_tag>()
    );
    FEDofMoFEMEntity_multiIndex::index<Composite_EntType_and_Space_mi_tag>::type::iterator dit,hi_dit;
    dit = data_dofs.lower_bound(boost::make_tuple(type,space));
    hi_dit = data_dofs.upper_bound(boost::make_tuple(type,space));
    for(;dit!=hi_dit;dit++) {
      ApproximationOrder ent_order = dit->get_max_order();
      int side_number = dit->side_number_ptr->side_number;
      if(side_number < 0) {
        SETERRQ(PETSC_COMM_SELF,MOFEM_DATA_INCONSISTENCY,"data inconsistency");
      }
      data[side_number].getDataOrder() = data[side_number].getDataOrder() > ent_order ? data[side_number].getDataOrder() : ent_order;
      if(dit->side_number_ptr->brother_side_number!=-1) {
        if(data.size() < (unsigned int)dit->side_number_ptr->brother_side_number) {
          SETERRQ(PETSC_COMM_SELF,MOFEM_DATA_INCONSISTENCY,"data inconsistency");
        }
        data[dit->side_number_ptr->brother_side_number].getDataOrder() = data[side_number].getDataOrder();
      }
    }
  } catch (exception& ex) {
    ostringstream ss;
    ss << "thorw in method: " << ex.what() << " at line " << __LINE__ << " in file " << __FILE__;
    SETERRQ(PETSC_COMM_SELF,MOFEM_STD_EXCEPTION_THROW,ss.str().c_str());
  }
  PetscFunctionReturn(0);
}

PetscErrorCode ForcesAndSurcesCore::getEdgesDataOrder(DataForcesAndSurcesCore &data,const FieldSpace space) {
  PetscFunctionBegin;
  ierr = getDataOrder(MBEDGE,space,data.dataOnEntities[MBEDGE]); CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

PetscErrorCode ForcesAndSurcesCore::getTrisDataOrder(DataForcesAndSurcesCore &data,const FieldSpace space) {
  PetscFunctionBegin;
  ierr = getDataOrder(MBTRI,space,data.dataOnEntities[MBTRI]); CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

PetscErrorCode ForcesAndSurcesCore::getQuadDataOrder(DataForcesAndSurcesCore &data,const FieldSpace space) {
  PetscFunctionBegin;
  ierr = getDataOrder(MBQUAD,space,data.dataOnEntities[MBQUAD]); CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

PetscErrorCode ForcesAndSurcesCore::getTetDataOrder(DataForcesAndSurcesCore &data,const FieldSpace space) {
  PetscFunctionBegin;
  ierr = getDataOrder(MBTET,space,data.dataOnEntities[MBTET]); CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

PetscErrorCode ForcesAndSurcesCore::getPrismDataOrder(DataForcesAndSurcesCore &data,const FieldSpace space) {
  PetscFunctionBegin;
  ierr = getDataOrder(MBPRISM,space,data.dataOnEntities[MBPRISM]); CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

PetscErrorCode ForcesAndSurcesCore::getDataOrderSpaceAndBase(
  const string &field_name,const EntityType type,boost::ptr_vector<DataForcesAndSurcesCore::EntData> &data
) {
  PetscFunctionBegin;

  SideNumber_multiIndex& side_table = const_cast<SideNumber_multiIndex&>(fePtr->get_side_number_table());
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
  }

  FEDofMoFEMEntity_multiIndex::index<Composite_Name_And_Type_mi_tag>::type &data_dofs =
  const_cast<FEDofMoFEMEntity_multiIndex::index<Composite_Name_And_Type_mi_tag>::type&>(
    fePtr->get_data_dofs().get<Composite_Name_And_Type_mi_tag>()
  );

  FEDofMoFEMEntity_multiIndex::index<Composite_Name_And_Type_mi_tag>::type::iterator dit,hi_dit;
  dit = data_dofs.lower_bound(boost::make_tuple(field_name,type));
  hi_dit = data_dofs.upper_bound(boost::make_tuple(field_name,type));

  for(;dit!=hi_dit;dit++) {

    // cerr << ApproximationBaseNames[dit->get_approx_base()] << endl;

    int side_number = dit->side_number_ptr->side_number;
    if(data[side_number].getDataOrder()) continue;

    ApproximationOrder ent_order = dit->get_max_order();
    if(side_number < 0) {
      SETERRQ(PETSC_COMM_SELF,MOFEM_DATA_INCONSISTENCY,"data inconsistency");
    }
    data[side_number].getBase() = dit->get_approx_base();
    data[side_number].getSpace() = dit->get_space();
    data[side_number].getDataOrder() = data[side_number].getDataOrder() > ent_order ? data[side_number].getDataOrder() : ent_order;
    if(dit->side_number_ptr->brother_side_number!=-1) {
      if(data.size() < (unsigned int)dit->side_number_ptr->brother_side_number) {
        SETERRQ(PETSC_COMM_SELF,MOFEM_DATA_INCONSISTENCY,"data inconsistency");
      }
      data[dit->side_number_ptr->brother_side_number].getBase() = data[side_number].getBase();
      data[dit->side_number_ptr->brother_side_number].getSpace() = data[side_number].getSpace();
      data[dit->side_number_ptr->brother_side_number].getDataOrder() = data[side_number].getDataOrder();
    }

  }

  PetscFunctionReturn(0);
}

PetscErrorCode ForcesAndSurcesCore::getEdgesDataOrderSpaceAndBase(DataForcesAndSurcesCore &data,const string &field_name) {
  PetscFunctionBegin;
  ierr = getDataOrderSpaceAndBase(field_name,MBEDGE,data.dataOnEntities[MBEDGE]); CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

PetscErrorCode ForcesAndSurcesCore::getTrisDataOrderSpaceAndBase(DataForcesAndSurcesCore &data,const string &field_name) {
  PetscFunctionBegin;
  ierr = getDataOrderSpaceAndBase(field_name,MBTRI,data.dataOnEntities[MBTRI]); CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

PetscErrorCode ForcesAndSurcesCore::getQuadDataOrderSpaceAndBase(DataForcesAndSurcesCore &data,const string &field_name) {
  PetscFunctionBegin;
  ierr = getDataOrderSpaceAndBase(field_name,MBQUAD,data.dataOnEntities[MBQUAD]); CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

PetscErrorCode ForcesAndSurcesCore::getTetDataOrderSpaceAndBase(DataForcesAndSurcesCore &data,const string &field_name) {
  PetscFunctionBegin;
  ierr = getDataOrderSpaceAndBase(field_name,MBTET,data.dataOnEntities[MBTET]); CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

PetscErrorCode ForcesAndSurcesCore::getPrismDataOrderSpaceAndBase(DataForcesAndSurcesCore &data,const string &field_name) {
  PetscFunctionBegin;
  ierr = getDataOrderSpaceAndBase(field_name,MBPRISM,data.dataOnEntities[MBPRISM]); CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

// ** Indices **

PetscErrorCode ForcesAndSurcesCore::getNodesIndices(
  const string &field_name,FENumeredDofMoFEMEntity_multiIndex &dofs,VectorInt &nodes_indices,VectorInt &local_nodes_indices
) {
  PetscFunctionBegin;
  FENumeredDofMoFEMEntity_multiIndex::index<Composite_Name_And_Type_mi_tag>::type::iterator dit,hi_dit,it;
  dit = dofs.get<Composite_Name_And_Type_mi_tag>().lower_bound(boost::make_tuple(field_name,MBVERTEX));
  hi_dit = dofs.get<Composite_Name_And_Type_mi_tag>().upper_bound(boost::make_tuple(field_name,MBVERTEX));

  int num_nodes;
  ierr = getNumberOfNodes(num_nodes); CHKERRQ(ierr);
  int max_nb_dofs = 0;
  if(dit!=hi_dit) {
    max_nb_dofs = dit->get_nb_of_coeffs()*num_nodes;
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
    int idx = dit->get_petsc_gloabl_dof_idx();
    int local_idx = dit->get_petsc_local_dof_idx();
    int side_number = dit->side_number_ptr->side_number;
    int pos = side_number*dit->get_nb_of_coeffs()+dit->get_dof_coeff_idx();
    nodes_indices[pos] = idx;
    local_nodes_indices[pos] = local_idx;
    int  brother_side_number = dit->side_number_ptr->brother_side_number;
    if(brother_side_number!=-1) {
      if(nodes_indices.size()<(unsigned int)(brother_side_number*dit->get_nb_of_coeffs()+dit->get_nb_of_coeffs())) {
        nodes_indices.resize(brother_side_number*dit->get_nb_of_coeffs()+dit->get_nb_of_coeffs());
      }
      int elem_idx = brother_side_number*dit->get_nb_of_coeffs()+dit->get_dof_coeff_idx();
      nodes_indices[elem_idx] = idx;
      local_nodes_indices[elem_idx] = local_idx;
    }
  }
  PetscFunctionReturn(0);
}

PetscErrorCode ForcesAndSurcesCore::getRowNodesIndices(DataForcesAndSurcesCore &data,const string &field_name) {
  PetscFunctionBegin;
  ierr = getNodesIndices(
    field_name,const_cast<FENumeredDofMoFEMEntity_multiIndex&>(fePtr->get_rows_dofs()),
    data.dataOnEntities[MBVERTEX][0].getIndices(),data.dataOnEntities[MBVERTEX][0].getLocalIndices()
  ); CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

PetscErrorCode ForcesAndSurcesCore::getColNodesIndices(DataForcesAndSurcesCore &data,const string &field_name) {
  PetscFunctionBegin;
  ierr = getNodesIndices(
    field_name,const_cast<FENumeredDofMoFEMEntity_multiIndex&>(fePtr->get_cols_dofs()),
    data.dataOnEntities[MBVERTEX][0].getIndices(),data.dataOnEntities[MBVERTEX][0].getLocalIndices()
  ); CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

PetscErrorCode ForcesAndSurcesCore::getTypeIndices(
  const string &field_name,FENumeredDofMoFEMEntity_multiIndex &dofs,EntityType type,int side_number,VectorInt &indices,VectorInt &local_indices
) {
  PetscFunctionBegin;
  FENumeredDofMoFEMEntity_multiIndex::index<Composite_Name_Type_And_Side_Number_mi_tag>::type::iterator dit,hi_dit;
  dit = dofs.get<Composite_Name_Type_And_Side_Number_mi_tag>().lower_bound(boost::make_tuple(field_name,type,side_number));
  hi_dit = dofs.get<Composite_Name_Type_And_Side_Number_mi_tag>().upper_bound(boost::make_tuple(field_name,type,side_number));
  if(dit!=hi_dit) {
    indices.resize(dit->get_nb_dofs_on_ent(),false);
    local_indices.resize(dit->get_nb_dofs_on_ent(),false);
    for(;dit!=hi_dit;dit++) {
      int idx = dit->get_petsc_gloabl_dof_idx();
      int elemem_idx = dit->get_EntDofIdx();
      indices[elemem_idx] = idx;
      int local_idx = dit->get_petsc_local_dof_idx();
      local_indices[elemem_idx] = local_idx;
    }
  } else {
    indices.resize(0);
    local_indices.resize(0);
  }
  PetscFunctionReturn(0);
}

PetscErrorCode ForcesAndSurcesCore::getTypeIndices(
  const string &field_name,FENumeredDofMoFEMEntity_multiIndex &dofs,EntityType type,
  boost::ptr_vector<DataForcesAndSurcesCore::EntData> &data
) {
  PetscFunctionBegin;
  SideNumber_multiIndex& side_table = const_cast<SideNumber_multiIndex&>(fePtr->get_side_number_table());
  SideNumber_multiIndex::nth_index<2>::type::iterator siit = side_table.get<2>().lower_bound(type);
  SideNumber_multiIndex::nth_index<2>::type::iterator hi_siit = side_table.get<2>().upper_bound(type);
  for(;siit!=hi_siit;siit++) {
    ierr = getTypeIndices(
      field_name,dofs,type,siit->side_number,data[siit->side_number].getIndices(),data[siit->side_number].getLocalIndices()
    ); CHKERRQ(ierr);
    if(siit->brother_side_number!=-1) {
      ierr = getTypeIndices(
        field_name,dofs,type,siit->side_number,data[siit->brother_side_number].getIndices(),data[siit->brother_side_number].getLocalIndices()
      ); CHKERRQ(ierr);
    }
  }
  PetscFunctionReturn(0);
}

PetscErrorCode ForcesAndSurcesCore::getEdgesRowIndices(DataForcesAndSurcesCore &data,const string &field_name) {
  PetscFunctionBegin;
  ierr = getTypeIndices(
    field_name,const_cast<FENumeredDofMoFEMEntity_multiIndex&>(fePtr->get_rows_dofs()),MBEDGE,data.dataOnEntities[MBEDGE]
  ); CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

PetscErrorCode ForcesAndSurcesCore::getEdgesColIndices(DataForcesAndSurcesCore &data,const string &field_name) {
  PetscFunctionBegin;
  ierr = getTypeIndices(
    field_name,const_cast<FENumeredDofMoFEMEntity_multiIndex&>(fePtr->get_cols_dofs()),MBEDGE,data.dataOnEntities[MBEDGE]
  ); CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

PetscErrorCode ForcesAndSurcesCore::getTrisRowIndices(
  DataForcesAndSurcesCore &data,const string &field_name
) {
  PetscFunctionBegin;
  ierr = getTypeIndices(
    field_name,const_cast<FENumeredDofMoFEMEntity_multiIndex&>(fePtr->get_rows_dofs()),MBTRI,data.dataOnEntities[MBTRI]
  ); CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

PetscErrorCode ForcesAndSurcesCore::getTrisColIndices(DataForcesAndSurcesCore &data,const string &field_name) {
  PetscFunctionBegin;
  ierr = getTypeIndices(
    field_name,const_cast<FENumeredDofMoFEMEntity_multiIndex&>(fePtr->get_cols_dofs()),MBTRI,data.dataOnEntities[MBTRI]
  ); CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

PetscErrorCode ForcesAndSurcesCore::getTetsRowIndices(DataForcesAndSurcesCore &data,const string &field_name) {
  PetscFunctionBegin;
  if(data.dataOnEntities[MBTET].size() == 0) {
    SETERRQ(PETSC_COMM_SELF,MOFEM_DATA_INCONSISTENCY,"data inconsistency");
  }
  ierr = getTypeIndices(
    field_name,const_cast<FENumeredDofMoFEMEntity_multiIndex&>(fePtr->get_rows_dofs()),MBTET,0,
    data.dataOnEntities[MBTET][0].getIndices(),data.dataOnEntities[MBTET][0].getLocalIndices()
  ); CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

PetscErrorCode ForcesAndSurcesCore::getTetsColIndices(DataForcesAndSurcesCore &data,const string &field_name) {
  PetscFunctionBegin;
  if(data.dataOnEntities[MBTET].size() == 0) {
    SETERRQ(PETSC_COMM_SELF,MOFEM_DATA_INCONSISTENCY,"data inconsistency");
  }
  ierr = getTypeIndices(
    field_name,const_cast<FENumeredDofMoFEMEntity_multiIndex&>(
      fePtr->get_cols_dofs()),MBTET,0,data.dataOnEntities[MBTET][0].getIndices(),data.dataOnEntities[MBTET][0].getLocalIndices()
  ); CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

PetscErrorCode ForcesAndSurcesCore::getQuadRowIndices(
  DataForcesAndSurcesCore &data,const string &field_name
) {
  PetscFunctionBegin;
  ierr = getTypeIndices(
    field_name,const_cast<FENumeredDofMoFEMEntity_multiIndex&>(fePtr->get_rows_dofs()),MBQUAD,data.dataOnEntities[MBQUAD]
  ); CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

PetscErrorCode ForcesAndSurcesCore::getQuadColIndices(DataForcesAndSurcesCore &data,const string &field_name) {
  PetscFunctionBegin;
  ierr = getTypeIndices(
    field_name,const_cast<FENumeredDofMoFEMEntity_multiIndex&>(fePtr->get_cols_dofs()),MBQUAD,data.dataOnEntities[MBQUAD]
  ); CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

PetscErrorCode ForcesAndSurcesCore::getPrismRowIndices(
  DataForcesAndSurcesCore &data,const string &field_name
) {
  PetscFunctionBegin;
  ierr = getTypeIndices(
    field_name,const_cast<FENumeredDofMoFEMEntity_multiIndex&>(fePtr->get_rows_dofs()),MBPRISM,data.dataOnEntities[MBPRISM]
  ); CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

PetscErrorCode ForcesAndSurcesCore::getPrismColIndices(DataForcesAndSurcesCore &data,const string &field_name) {
  PetscFunctionBegin;
  ierr = getTypeIndices(
    field_name,const_cast<FENumeredDofMoFEMEntity_multiIndex&>(fePtr->get_cols_dofs()),MBPRISM,data.dataOnEntities[MBPRISM]
  ); CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

PetscErrorCode ForcesAndSurcesCore::getNoFieldIndices(
  const string &field_name,FENumeredDofMoFEMEntity_multiIndex &dofs,VectorInt &indices
) {
  PetscFunctionBegin;
  FENumeredDofMoFEMEntity_multiIndex::index<FieldName_mi_tag>::type::iterator dit,hi_dit;
  dit = dofs.get<FieldName_mi_tag>().lower_bound(field_name);
  hi_dit = dofs.get<FieldName_mi_tag>().upper_bound(field_name);
  indices.resize(distance(dit,hi_dit));
  for(;dit!=hi_dit;dit++) {
    int idx = dit->get_petsc_gloabl_dof_idx();
    indices[dit->get_dof_coeff_idx()] = idx;
  }
  PetscFunctionReturn(0);
}

PetscErrorCode ForcesAndSurcesCore::getNoFieldRowIndices(DataForcesAndSurcesCore &data,const string &field_name) {
  PetscFunctionBegin;
  //EntityType fe_type = fePtr->get_ent_type();
  if(data.dataOnEntities[MBENTITYSET].size() == 0) {
    SETERRQ(PETSC_COMM_SELF,MOFEM_DATA_INCONSISTENCY,"data inconsistency");
  }
  ierr = getNoFieldIndices(
    field_name,const_cast<FENumeredDofMoFEMEntity_multiIndex&>(fePtr->get_rows_dofs()),data.dataOnEntities[MBENTITYSET][0].getIndices()
  ); CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

PetscErrorCode ForcesAndSurcesCore::getNoFieldColIndices(DataForcesAndSurcesCore &data,const string &field_name) {
  PetscFunctionBegin;
  if(data.dataOnEntities[MBENTITYSET].size() == 0) {
    SETERRQ(PETSC_COMM_SELF,MOFEM_DATA_INCONSISTENCY,"data inconsistency");
  }
  ierr = getNoFieldIndices(
    field_name,const_cast<FENumeredDofMoFEMEntity_multiIndex&>(fePtr->get_cols_dofs()),data.dataOnEntities[MBENTITYSET][0].getIndices()
  ); CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

// ** Indices from problem **

PetscErrorCode ForcesAndSurcesCore::getProblemNodesIndices(const string &field_name,
  const NumeredDofMoFEMEntity_multiIndex &dofs,
  VectorInt &nodes_indices
) const {
  PetscFunctionBegin;

  nodes_indices.resize(0);

  SideNumber_multiIndex& side_table = const_cast<SideNumber_multiIndex&>(fePtr->get_side_number_table());
  SideNumber_multiIndex::nth_index<1>::type::iterator siit = side_table.get<1>().lower_bound(boost::make_tuple(MBVERTEX,0));
  SideNumber_multiIndex::nth_index<1>::type::iterator hi_siit = side_table.get<1>().lower_bound(boost::make_tuple(MBVERTEX,10000));

  int nn = 0;
  for(;siit!=hi_siit;siit++,nn++) {

    if(siit->side_number == -1) continue;

    const EntityHandle ent = siit->ent;
    NumeredDofMoFEMEntity_multiIndex::index<Composite_Name_And_Ent_And_EndDofIdx_mi_tag>::type::iterator dit,hi_dit;
    dit = dofs.get<Composite_Name_And_Ent_And_EndDofIdx_mi_tag>().lower_bound(boost::make_tuple(field_name,ent,0));
    hi_dit = dofs.get<Composite_Name_And_Ent_And_EndDofIdx_mi_tag>().upper_bound(boost::make_tuple(field_name,ent,10000));  /// very large number

    if(dit!=hi_dit) {

      if(!nn) {
        nodes_indices.resize(dit->get_nb_of_coeffs()*distance(siit,hi_siit));
      }
      for(;dit!=hi_dit;dit++) {
        nodes_indices[siit->side_number*dit->get_nb_of_coeffs()+dit->get_dof_coeff_idx()] = dit->get_petsc_gloabl_dof_idx();
      }

    }

  }

  PetscFunctionReturn(0);

}

PetscErrorCode ForcesAndSurcesCore::getProblemTypeIndices(
  const string &field_name,const NumeredDofMoFEMEntity_multiIndex &dofs,
  EntityType type,int side_number,VectorInt &indices
) const {
  PetscFunctionBegin;

  indices.resize(0);

  SideNumber_multiIndex& side_table = const_cast<SideNumber_multiIndex&>(fePtr->get_side_number_table());
  SideNumber_multiIndex::nth_index<1>::type::iterator siit = side_table.get<1>().lower_bound(boost::make_tuple(type,side_number));
  SideNumber_multiIndex::nth_index<1>::type::iterator hi_siit = side_table.get<1>().upper_bound(boost::make_tuple(type,side_number));

  for(;siit!=hi_siit;siit++) {

    if(siit->side_number == -1) continue;

    const EntityHandle ent = siit->ent;
    NumeredDofMoFEMEntity_multiIndex::index<Composite_Name_And_Ent_And_EndDofIdx_mi_tag>::type::iterator dit,hi_dit;
    dit = dofs.get<Composite_Name_And_Ent_And_EndDofIdx_mi_tag>().lower_bound(boost::make_tuple(field_name,ent,0));
    hi_dit = dofs.get<Composite_Name_And_Ent_And_EndDofIdx_mi_tag>().upper_bound(boost::make_tuple(field_name,ent,10000));  /// very large number

    indices.resize(distance(dit,hi_dit));
    for(;dit!=hi_dit;dit++) {

      indices[dit->get_EntDofIdx()] = dit->get_petsc_gloabl_dof_idx();

    }

  }

  PetscFunctionReturn(0);
}

PetscErrorCode ForcesAndSurcesCore::getProblemNodesRowIndices(const string &field_name,VectorInt &nodes_indices) const {
  PetscFunctionBegin;
  PetscErrorCode ierr;
  ierr = getProblemNodesIndices(field_name,problemPtr->numered_dofs_rows,nodes_indices); CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

PetscErrorCode ForcesAndSurcesCore::getProblemTypeRowIndices(const string &field_name,EntityType type,int side_number,VectorInt &indices) const {
  PetscFunctionBegin;
  PetscErrorCode ierr;
  ierr = getProblemTypeIndices(field_name,problemPtr->numered_dofs_rows,type,side_number,indices); CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

PetscErrorCode ForcesAndSurcesCore::getProblemNodesColIndices(const string &field_name,VectorInt &nodes_indices) const {
  PetscFunctionBegin;
  PetscErrorCode ierr;
  ierr = getProblemNodesIndices(field_name,problemPtr->numered_dofs_cols,nodes_indices); CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

PetscErrorCode ForcesAndSurcesCore::getProblemTypeColIndices(const string &field_name,EntityType type,int side_number,VectorInt &indices) const {
  PetscFunctionBegin;
  PetscErrorCode ierr;
  ierr = getProblemTypeIndices(field_name,problemPtr->numered_dofs_cols,type,side_number,indices); CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

// ** Data **

PetscErrorCode ForcesAndSurcesCore::getNodesFieldData(
  const string &field_name,
  FEDofMoFEMEntity_multiIndex &dofs,
  VectorDouble &nodes_data,
  VectorDofs &nodes_dofs,
  FieldSpace &space,
  FieldApproximationBase &base
) {
  PetscFunctionBegin;
  try {
    FEDofMoFEMEntity_multiIndex::index<Composite_Name_And_Type_mi_tag>::type::iterator dit,hi_dit,it;
    dit = dofs.get<Composite_Name_And_Type_mi_tag>().lower_bound(boost::make_tuple(field_name,MBVERTEX));
    hi_dit = dofs.get<Composite_Name_And_Type_mi_tag>().upper_bound(boost::make_tuple(field_name,MBVERTEX));

    int num_nodes;
    ierr = getNumberOfNodes(num_nodes); CHKERRQ(ierr);
    int max_nb_dofs = 0;
    if(dit!=hi_dit) {
      max_nb_dofs = dit->get_nb_of_coeffs()*num_nodes;
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
      space = dit->get_space();
      base = dit->get_approx_base();
    }

    for(;dit!=hi_dit;dit++) {
      FieldData val = dit->get_FieldData();
      int side_number = dit->side_number_ptr->side_number;
      if(side_number == -1) {
        SETERRQ(PETSC_COMM_SELF,MOFEM_DATA_INCONSISTENCY,"data inconsistency");
      }
      int pos = side_number*dit->get_nb_of_coeffs()+dit->get_dof_coeff_idx();
      nodes_data[pos] = val;
      nodes_dofs[pos] = &*dit;
      int  brother_side_number = dit->side_number_ptr->brother_side_number;
      if(brother_side_number!=-1) {
        if(nodes_data.size()<(unsigned int)(brother_side_number*dit->get_nb_of_coeffs()+dit->get_nb_of_coeffs())) {
          nodes_data.resize(brother_side_number*dit->get_nb_of_coeffs()+dit->get_nb_of_coeffs());
          nodes_dofs.resize(brother_side_number*dit->get_nb_of_coeffs()+dit->get_nb_of_coeffs());
        }
        int brother_pos = brother_side_number*dit->get_nb_of_coeffs()+dit->get_dof_coeff_idx();
        nodes_data[brother_pos] = val;
        nodes_dofs[brother_pos] = &*dit;
      }
    }

  } catch (exception& ex) {
    ostringstream ss;
    ss << "thorw in method: " << ex.what() << " at line " << __LINE__ << " in file " << __FILE__;
    SETERRQ(PETSC_COMM_SELF,MOFEM_STD_EXCEPTION_THROW,ss.str().c_str());
  }
  PetscFunctionReturn(0);
}

PetscErrorCode ForcesAndSurcesCore::getNodesFieldData(DataForcesAndSurcesCore &data,const string &field_name) {
  PetscFunctionBegin;
  ierr = getNodesFieldData(
    field_name,
    const_cast<FEDofMoFEMEntity_multiIndex&>(fePtr->get_data_dofs()),
    data.dataOnEntities[MBVERTEX][0].getFieldData(),
    data.dataOnEntities[MBVERTEX][0].getFieldDofs(),
    data.dataOnEntities[MBVERTEX][0].getSpace(),
    data.dataOnEntities[MBVERTEX][0].getBase()
  ); CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

PetscErrorCode ForcesAndSurcesCore::getTypeFieldData(
  const string &field_name,
  FEDofMoFEMEntity_multiIndex &dofs,
  EntityType type,
  int side_number,
  VectorDouble &ent_field_data,
  VectorDofs &ent_field_dofs
) {
  PetscFunctionBegin;
  FEDofMoFEMEntity_multiIndex::index<Composite_Name_Type_And_Side_Number_mi_tag>::type::iterator dit,hi_dit;
  dit = dofs.get<Composite_Name_Type_And_Side_Number_mi_tag>().lower_bound(boost::make_tuple(field_name,type,side_number));
  hi_dit = dofs.get<Composite_Name_Type_And_Side_Number_mi_tag>().upper_bound(boost::make_tuple(field_name,type,side_number));
  if(dit!=hi_dit) {
    ent_field_data.resize(dit->get_nb_dofs_on_ent(),false);
    ent_field_dofs.resize(dit->get_nb_dofs_on_ent(),false);
    for(;dit!=hi_dit;dit++) {
      FieldData val = dit->get_FieldData();
      int idx = dit->get_EntDofIdx();
      ent_field_data[idx] = val;
      ent_field_dofs[idx] = &*dit;
    }
  } else {
    ent_field_data.resize(0);
    ent_field_dofs.resize(0);
  }
  PetscFunctionReturn(0);
}

PetscErrorCode ForcesAndSurcesCore::getTypeFieldData(
  const string &field_name,
  FEDofMoFEMEntity_multiIndex &dofs,
  EntityType type,
  boost::ptr_vector<DataForcesAndSurcesCore::EntData> &data
) {
  PetscFunctionBegin;
  SideNumber_multiIndex& side_table = const_cast<SideNumber_multiIndex&>(fePtr->get_side_number_table());
  if(data.size() < side_table.get<2>().count(type)) {
    SETERRQ(PETSC_COMM_SELF,MOFEM_DATA_INCONSISTENCY,"data inconsistency");
  }
  SideNumber_multiIndex::nth_index<2>::type::iterator siit = side_table.get<2>().lower_bound(type);
  SideNumber_multiIndex::nth_index<2>::type::iterator hi_siit = side_table.get<2>().upper_bound(type);
  for(;siit!=hi_siit;siit++) {
    ierr = getTypeFieldData(
      field_name,dofs,type,
      siit->side_number,
      data[siit->side_number].getFieldData(),
      data[siit->side_number].getFieldDofs()
    ); CHKERRQ(ierr);
    if(siit->brother_side_number!=-1) {
      if(data.size() < (unsigned int)siit->brother_side_number) {
        SETERRQ(PETSC_COMM_SELF,MOFEM_DATA_INCONSISTENCY,"data inconsistency");
      }
      ierr = getTypeFieldData(
        field_name,dofs,type,siit->side_number,
        data[siit->brother_side_number].getFieldData(),
        data[siit->brother_side_number].getFieldDofs()
      ); CHKERRQ(ierr);
    }
  }
  PetscFunctionReturn(0);
}

PetscErrorCode ForcesAndSurcesCore::getNoFieldFieldData(
  const string &field_name,
  FEDofMoFEMEntity_multiIndex &dofs,
  VectorDouble &ent_field_data,
  VectorDofs &ent_field_dofs
) {
  PetscFunctionBegin;
  FEDofMoFEMEntity_multiIndex::index<FieldName_mi_tag>::type::iterator dit,hi_dit;
  dit = dofs.get<FieldName_mi_tag>().lower_bound(field_name);
  hi_dit = dofs.get<FieldName_mi_tag>().upper_bound(field_name);
  int size = distance(dit,hi_dit);
  ent_field_data.resize(size,false);
  ent_field_dofs.resize(size,false);
  for(;dit!=hi_dit;dit++) {
    int idx = dit->get_dof_coeff_idx();
    ent_field_data[idx] = dit->get_FieldData();
    ent_field_dofs[idx] = &*dit;
  }
  PetscFunctionReturn(0);
}

PetscErrorCode ForcesAndSurcesCore::getNoFieldFieldData(
  DataForcesAndSurcesCore &data,const string &field_name
) {
  PetscFunctionBegin;
  if(data.dataOnEntities[MBENTITYSET].size() == 0) {
    SETERRQ(PETSC_COMM_SELF,MOFEM_DATA_INCONSISTENCY,"data inconsistency");
  }
  ierr = getNoFieldFieldData(
    field_name,const_cast<FEDofMoFEMEntity_multiIndex&>(fePtr->get_data_dofs()),
    data.dataOnEntities[MBENTITYSET][0].getFieldData(),
    data.dataOnEntities[MBENTITYSET][0].getFieldDofs()
  ); CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

PetscErrorCode ForcesAndSurcesCore::getEdgesFieldData(DataForcesAndSurcesCore &data,const string &field_name) {
  PetscFunctionBegin;
  ierr = getTypeFieldData(
    field_name,const_cast<FEDofMoFEMEntity_multiIndex&>(fePtr->get_data_dofs()),MBEDGE,data.dataOnEntities[MBEDGE]
  ); CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

PetscErrorCode ForcesAndSurcesCore::getTrisFieldData(
  DataForcesAndSurcesCore &data,const string &field_name
) {
  PetscFunctionBegin;
  ierr = getTypeFieldData(
    field_name,const_cast<FEDofMoFEMEntity_multiIndex&>(fePtr->get_data_dofs()),MBTRI,data.dataOnEntities[MBTRI]
  ); CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

PetscErrorCode ForcesAndSurcesCore::getQuadFieldData(
  DataForcesAndSurcesCore &data,const string &field_name
) {
  PetscFunctionBegin;
  ierr = getTypeFieldData(
    field_name,const_cast<FEDofMoFEMEntity_multiIndex&>(fePtr->get_data_dofs()),MBQUAD,data.dataOnEntities[MBQUAD]
  ); CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

PetscErrorCode ForcesAndSurcesCore::getTetsFieldData(DataForcesAndSurcesCore &data,const string &field_name) {
  PetscFunctionBegin;
  if(data.dataOnEntities[MBTET].size() == 0) {
    SETERRQ(PETSC_COMM_SELF,MOFEM_DATA_INCONSISTENCY,"data inconsistency");
  }
  ierr = getTypeFieldData(
    field_name,const_cast<FEDofMoFEMEntity_multiIndex&>(fePtr->get_data_dofs()),
    MBTET,
    0,
    data.dataOnEntities[MBTET][0].getFieldData(),
    data.dataOnEntities[MBTET][0].getFieldDofs()
  ); CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

PetscErrorCode ForcesAndSurcesCore::getPrismFieldData(DataForcesAndSurcesCore &data,const string &field_name) {
  PetscFunctionBegin;
  if(data.dataOnEntities[MBPRISM].size() == 0) {
    SETERRQ(PETSC_COMM_SELF,MOFEM_DATA_INCONSISTENCY,"data inconsistency");
  }
  ierr = getTypeFieldData(
    field_name,
    const_cast<FEDofMoFEMEntity_multiIndex&>(fePtr->get_data_dofs()),
    MBPRISM,
    0,
    data.dataOnEntities[MBPRISM][0].getFieldData(),
    data.dataOnEntities[MBPRISM][0].getFieldDofs()
  ); CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

// ** Face **

PetscErrorCode ForcesAndSurcesCore::getFaceTriNodes(DataForcesAndSurcesCore &data) {
  PetscFunctionBegin;
  //PetscAttachDebugger();
  data.facesNodes.resize(4,3,false);
  SideNumber_multiIndex& side_table = const_cast<SideNumber_multiIndex&>(fePtr->get_side_number_table());
  SideNumber_multiIndex::nth_index<1>::type::iterator siit = side_table.get<1>().lower_bound(boost::make_tuple(MBTRI,0));
  SideNumber_multiIndex::nth_index<1>::type::iterator hi_siit = side_table.get<1>().upper_bound(boost::make_tuple(MBTRI,4));
  const int cannonical_face_sense_p1[4][3] = { {0,1,3}, {1,2,3}, {0,3,2}/**/, {0,2,1}/**/ }; //secon index is offset (positive sense)
  const int cannonical_face_sense_m1[4][3] = { {0,3,1}, {1,3,2}, {0,2,3}, {0,1,2} }; //second index is offset (negative sense
  for(;siit!=hi_siit;siit++) {
    const SideNumber* side = &*siit;
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
      EntityHandle ent = fePtr->get_ent();
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

PetscErrorCode ForcesAndSurcesCore::getSpacesAndBaseOnEntities(DataForcesAndSurcesCore &data) {
  PetscFunctionBegin;
  try {
    for(_IT_GET_FEDATA_DOFS_FOR_LOOP_(this,dof)) {
      // cerr << *dof << endl;
      // cerr << dof->get_space() << " " << data.sPace.size() << endl;
      // cerr << dof->get_approx_base() << " " << data.bAse.size() << endl;
      data.sPace.set(dof->get_space());
      data.bAse.set(dof->get_approx_base());
      data.spacesOnEntities[dof->get_ent_type()].set(dof->get_space());
      data.basesOnEntities[dof->get_ent_type()].set(dof->get_approx_base());
      // cerr << "approx base " << ApproximationBaseNames[dof->get_approx_base()] << " " << data.basesOnEntities[dof->get_ent_type()] << endl;
    }
  } catch (exception& ex) {
    ostringstream ss;
    ss << "thorw in method: " << ex.what() << " at line " << __LINE__ << " in file " << __FILE__;
    SETERRQ(PETSC_COMM_SELF,MOFEM_STD_EXCEPTION_THROW,ss.str().c_str());
  }
  PetscFunctionReturn(0);
}

// ** Shape Functions **

PetscErrorCode ForcesAndSurcesCore::shapeTETFunctions_H1(
  DataForcesAndSurcesCore &data,
  const double *G_X,
  const double *G_Y,
  const double *G_Z,
  const int G_DIM,
  const FieldApproximationBase base,
  PetscErrorCode (*base_polynomials)(int p,double s,double *diff_s,double *L,double *diffL,const int dim)
) {
  PetscFunctionBegin;

  if(data.dataOnEntities[MBVERTEX][0].getN(base).size2()!=4) {
    SETERRQ1(
      PETSC_COMM_SELF,MOFEM_DATA_INCONSISTENCY,"Base functions on nodes not set for base %s",ApproximationBaseNames[base]
    );
  }
  if(data.dataOnEntities[MBVERTEX][0].getN(base).size1()!=G_DIM) {
    SETERRQ1(
      PETSC_COMM_SELF,
      MOFEM_DATA_INCONSISTENCY,
      "Base functions or nodes has wrong number of integration points for base %s",
      ApproximationBaseNames[base]
    );
  }

  data.dataOnEntities[MBVERTEX][0].getBase() = base;
  data.dataOnEntities[MBVERTEX][0].getDiffN(base).resize(4,3,false);
  ierr = ShapeDiffMBTET(&*data.dataOnEntities[MBVERTEX][0].getDiffN(base).data().begin()); CHKERRQ(ierr);

  int _sense_[6],_order_[6];
  if(data.spacesOnEntities[MBEDGE].test(H1)) {
    //edges
    if(data.dataOnEntities[MBEDGE].size()!=6) {
      SETERRQ(PETSC_COMM_SELF,MOFEM_DATA_INCONSISTENCY,"data inconsistency");
    }
    double *_H1edgeN_[6],*_diffH1edgeN_[6];
    for(int ee = 0;ee<6;ee++) {
      data.dataOnEntities[MBEDGE][ee].getBase() = base;
      if(data.dataOnEntities[MBEDGE][ee].getSense() == 0) {
        SETERRQ(PETSC_COMM_SELF,MOFEM_DATA_INCONSISTENCY,"data inconsistency");
      }
      _sense_[ee] = data.dataOnEntities[MBEDGE][ee].getSense();
      _order_[ee] = data.dataOnEntities[MBEDGE][ee].getDataOrder();
      int nb_dofs = NBEDGE_H1_AINSWORTH_COLE(data.dataOnEntities[MBEDGE][ee].getDataOrder());
      data.dataOnEntities[MBEDGE][ee].getN(base).resize(G_DIM,nb_dofs,false);
      data.dataOnEntities[MBEDGE][ee].getDiffN(base).resize(G_DIM,3*nb_dofs,false);
      _H1edgeN_[ee] = &*data.dataOnEntities[MBEDGE][ee].getN(base).data().begin();
      _diffH1edgeN_[ee] = &*data.dataOnEntities[MBEDGE][ee].getDiffN(base).data().begin();
    }
    ierr = H1_EdgeShapeFunctions_MBTET(
      _sense_,_order_,
      &*data.dataOnEntities[MBVERTEX][0].getN(base).data().begin(),
      &*data.dataOnEntities[MBVERTEX][0].getDiffN(base).data().begin(),
      _H1edgeN_,_diffH1edgeN_,G_DIM,base_polynomials
    ); CHKERRQ(ierr);
  }

  if(data.spacesOnEntities[MBTRI].test(H1)) {

    //faces
    if(data.dataOnEntities[MBTRI].size()!=4) {
      SETERRQ(PETSC_COMM_SELF,MOFEM_DATA_INCONSISTENCY,"data inconsistency");
    }
    double *_H1faceN_[4],*_diffH1faceN_[4];
    for(int ff = 0;ff<4;ff++) {
      data.dataOnEntities[MBTRI][ff].getBase() = base;
      if(data.dataOnEntities[MBTRI][ff].getSense() == 0) {
        SETERRQ(PETSC_COMM_SELF,MOFEM_DATA_INCONSISTENCY,"data inconsistency");
      }
      int nb_dofs = NBFACETRI_H1_AINSWORTH_COLE(data.dataOnEntities[MBTRI][ff].getDataOrder());
      _order_[ff] = data.dataOnEntities[MBTRI][ff].getDataOrder();
      data.dataOnEntities[MBTRI][ff].getN(base).resize(G_DIM,nb_dofs,false);
      data.dataOnEntities[MBTRI][ff].getDiffN(base).resize(G_DIM,3*nb_dofs,false);
      _H1faceN_[ff] = &*data.dataOnEntities[MBTRI][ff].getN(base).data().begin();
      _diffH1faceN_[ff] = &*data.dataOnEntities[MBTRI][ff].getDiffN(base).data().begin();
    }
    if(data.facesNodes.size1() != 4) {
      SETERRQ(PETSC_COMM_SELF,MOFEM_DATA_INCONSISTENCY,"data inconsistency");
    }
    if(data.facesNodes.size2() != 3) {
      SETERRQ(PETSC_COMM_SELF,MOFEM_DATA_INCONSISTENCY,"data inconsistency");
    }
    ierr = H1_FaceShapeFunctions_MBTET(
      &*data.facesNodes.data().begin(),
      _order_,
      &*data.dataOnEntities[MBVERTEX][0].getN(base).data().begin(),
      &*data.dataOnEntities[MBVERTEX][0].getDiffN(base).data().begin(),
      _H1faceN_,
      _diffH1faceN_,
      G_DIM,
      base_polynomials
    ); CHKERRQ(ierr);

  }

  if(data.spacesOnEntities[MBTET].test(H1)) {

    //volume
    data.dataOnEntities[MBTET][0].getBase() = base;
    int order = data.dataOnEntities[MBTET][0].getDataOrder();
    int nb_vol_dofs = NBVOLUMETET_H1_AINSWORTH_COLE(order);
    data.dataOnEntities[MBTET][0].getN(base).resize(G_DIM,nb_vol_dofs,false);
    data.dataOnEntities[MBTET][0].getDiffN(base).resize(G_DIM,3*nb_vol_dofs,false);
    ierr = H1_VolumeShapeFunctions_MBTET(
      data.dataOnEntities[MBTET][0].getDataOrder(),
      &*data.dataOnEntities[MBVERTEX][0].getN(base).data().begin(),
      &*data.dataOnEntities[MBVERTEX][0].getDiffN(base).data().begin(),
      &*data.dataOnEntities[MBTET][0].getN(base).data().begin(),
      &*data.dataOnEntities[MBTET][0].getDiffN(base).data().begin(),
      G_DIM,
      base_polynomials
    ); CHKERRQ(ierr);

  }

  PetscFunctionReturn(0);
}

PetscErrorCode ForcesAndSurcesCore::shapeTETFunctions_L2(
  DataForcesAndSurcesCore &data,
  const double *G_X,
  const double *G_Y,
  const double *G_Z,
  const int G_DIM,
  const FieldApproximationBase base,
  PetscErrorCode (*base_polynomials)(int p,double s,double *diff_s,double *L,double *diffL,const int dim)
) {
  PetscFunctionBegin;

  data.dataOnEntities[MBVERTEX][0].getN(base).resize(G_DIM,4,false);
  ierr = ShapeMBTET(&*data.dataOnEntities[MBVERTEX][0].getN(base).data().begin(),G_X,G_Y,G_Z,G_DIM); CHKERRQ(ierr);
  data.dataOnEntities[MBVERTEX][0].getDiffN(base).resize(4,3,false);
  ierr = ShapeDiffMBTET(&*data.dataOnEntities[MBVERTEX][0].getDiffN(base).data().begin()); CHKERRQ(ierr);

  data.dataOnEntities[MBTET][0].getN(base).resize(G_DIM,NBVOLUMETET_L2_AINSWORTH_COLE(data.dataOnEntities[MBTET][0].getDataOrder()),false);
  data.dataOnEntities[MBTET][0].getDiffN(base).resize(G_DIM,3*NBVOLUMETET_L2_AINSWORTH_COLE(data.dataOnEntities[MBTET][0].getDataOrder()),false);

  ierr = L2_ShapeFunctions_MBTET(
    data.dataOnEntities[MBTET][0].getDataOrder(),
    &*data.dataOnEntities[MBVERTEX][0].getN(base).data().begin(),
    &*data.dataOnEntities[MBVERTEX][0].getDiffN(base).data().begin(),
    &*data.dataOnEntities[MBTET][0].getN(base).data().begin(),
    &*data.dataOnEntities[MBTET][0].getDiffN(base).data().begin(),
    G_DIM
  ); CHKERRQ(ierr);

  PetscFunctionReturn(0);
}

PetscErrorCode ForcesAndSurcesCore::shapeTETFunctions_Hdiv(
  DataForcesAndSurcesCore &data,
  const double *G_X,
  const double *G_Y,
  const double *G_Z,
  const int G_DIM,
  const FieldApproximationBase base,
  PetscErrorCode (*base_polynomials)(int p,double s,double *diff_s,double *L,double *diffL,const int dim)
) {
  PetscFunctionBegin;

  //calculate shape function for tet, needed to construct shape functions for h_div

  if(data.dataOnEntities[MBVERTEX][0].getN(base).size1() != (unsigned int)G_DIM) {
    data.dataOnEntities[MBVERTEX][0].getN(base).resize(G_DIM,4,false);
    ierr = ShapeMBTET(&*data.dataOnEntities[MBVERTEX][0].getN(base).data().begin(),G_X,G_Y,G_Z,G_DIM); CHKERRQ(ierr);
  }
  //that is cheep to calate, no harm done if recalculated
  data.dataOnEntities[MBVERTEX][0].getDiffN(base).resize(4,3);
  ierr = ShapeDiffMBTET(&*data.dataOnEntities[MBVERTEX][0].getDiffN(base).data().begin()); CHKERRQ(ierr);

  //face shape functions

  double *phi_f_e[4][3];
  double *phi_f[4];
  double *diff_phi_f_e[4][3];
  double *diff_phi_f[4];

  N_face_edge.resize(4,3,false);
  N_face_bubble.resize(4,false);
  diffN_face_edge.resize(4,3,false);
  diffN_face_bubble.resize(4,false);

  int faces_order[4];
  for(int ff = 0;ff<4;ff++) {
    if(data.dataOnEntities[MBTRI][ff].getSense() == 0) {
      SETERRQ(PETSC_COMM_SELF,MOFEM_DATA_INCONSISTENCY,"data inconsistency");
    }
    faces_order[ff] = data.dataOnEntities[MBTRI][ff].getDataOrder();
    //three edges on face
    for(int ee = 0;ee<3;ee++) {
      N_face_edge(ff,ee).resize(G_DIM,3*NBFACETRI_EDGE_HDIV_AINSWORTH_COLE(faces_order[ff]),false);
      diffN_face_edge(ff,ee).resize(G_DIM,9*NBFACETRI_EDGE_HDIV_AINSWORTH_COLE(faces_order[ff]),false);
      phi_f_e[ff][ee] = &((N_face_edge(ff,ee))(0,0));
      diff_phi_f_e[ff][ee] = &((diffN_face_edge(ff,ee))(0,0));
    }
    N_face_bubble[ff].resize(G_DIM,3*NBFACETRI_FACE_HDIV_AINSWORTH_COLE(faces_order[ff]),false);
    diffN_face_bubble[ff].resize(G_DIM,9*NBFACETRI_FACE_HDIV_AINSWORTH_COLE(faces_order[ff]),false);
    phi_f[ff] = &*(N_face_bubble[ff].data().begin());
    diff_phi_f[ff] = &*(diffN_face_bubble[ff].data().begin());
  }

  ierr = Hdiv_EdgeFaceShapeFunctions_MBTET(
    &data.facesNodes(0,0),faces_order,
    &data.dataOnEntities[MBVERTEX][0].getN(base)(0,0),
    &data.dataOnEntities[MBVERTEX][0].getDiffN(base)(0,0),
    phi_f_e,diff_phi_f_e,G_DIM,
    base_polynomials
  ); CHKERRQ(ierr);

  ierr = Hdiv_FaceBubbleShapeFunctions_MBTET(
    &data.facesNodes(0,0),faces_order,
    &data.dataOnEntities[MBVERTEX][0].getN(base)(0,0),
    &data.dataOnEntities[MBVERTEX][0].getDiffN(base)(0,0),
    phi_f,diff_phi_f,G_DIM,
    base_polynomials
  ); CHKERRQ(ierr);

  //volume shape functions

  double *phi_v_e[6];
  double *phi_v_f[4];
  double *phi_v;
  double *diff_phi_v_e[6];
  double *diff_phi_v_f[4];
  double *diff_phi_v;

  int volume_order = data.dataOnEntities[MBTET][0].getDataOrder();
  double coords[] = { 0,0,0, 1,0,0, 0,1,0, 0,0,1 };

  N_volume_edge.resize(6,false);
  diffN_volume_edge.resize(6,false);
  for(int ee = 0;ee<6;ee++) {
    N_volume_edge[ee].resize(G_DIM,3*NBVOLUMETET_EDGE_HDIV_AINSWORTH_COLE(volume_order),false);
    diffN_volume_edge[ee].resize(G_DIM,9*NBVOLUMETET_EDGE_HDIV_AINSWORTH_COLE(volume_order),false);
    phi_v_e[ee] = &*(N_volume_edge[ee].data().begin());
    diff_phi_v_e[ee] = &*(diffN_volume_edge[ee].data().begin());
  }
  ierr = Hdiv_EdgeBasedVolumeShapeFunctions_MBTET(
    volume_order,coords,&data.dataOnEntities[MBVERTEX][0].getN(base)(0,0),
    &data.dataOnEntities[MBVERTEX][0].getDiffN(base)(0,0),
    phi_v_e,diff_phi_v_e,G_DIM,
    base_polynomials
  ); CHKERRQ(ierr);

  N_volume_face.resize(4,false);
  diffN_volume_face.resize(4,false);
  for(int ff = 0;ff<4;ff++) {
    N_volume_face[ff].resize(G_DIM,3*NBVOLUMETET_FACE_HDIV_AINSWORTH_COLE(volume_order),false);
    diffN_volume_face[ff].resize(G_DIM,9*NBVOLUMETET_FACE_HDIV_AINSWORTH_COLE(volume_order),false);
    phi_v_f[ff] = &*(N_volume_face[ff].data().begin());
    diff_phi_v_f[ff] = &*(diffN_volume_face[ff].data().begin());
  }
  ierr = Hdiv_FaceBasedVolumeShapeFunctions_MBTET(
    volume_order,coords,&data.dataOnEntities[MBVERTEX][0].getN(base)(0,0),
    &data.dataOnEntities[MBVERTEX][0].getDiffN(base)(0,0),
    phi_v_f,diff_phi_v_f,G_DIM,
    base_polynomials
  ); CHKERRQ(ierr);

  N_volume_bubble.resize(G_DIM,3*NBVOLUMETET_VOLUME_HDIV_AINSWORTH_COLE(volume_order),false);
  diffN_volume_bubble.resize(G_DIM,9*NBVOLUMETET_VOLUME_HDIV_AINSWORTH_COLE(volume_order),false);
  phi_v = &*(N_volume_bubble.data().begin());
  diff_phi_v = &*(diffN_volume_bubble.data().begin());
  ierr = Hdiv_VolumeBubbleShapeFunctions_MBTET(
    volume_order,coords,&data.dataOnEntities[MBVERTEX][0].getN(base)(0,0),
    &data.dataOnEntities[MBVERTEX][0].getDiffN(base)(0,0),
    phi_v,diff_phi_v,G_DIM,
    base_polynomials
  ); CHKERRQ(ierr);

  // Set shape functions into data strucrure Shape functions hast to be put
  // in arrays in order which guarantee hierarhical series of digrees of
  // freedom, i.e. in other words dofs form sub-entities has to be group
  // by order.

  //faces
  if(data.dataOnEntities[MBTRI].size()!=4) {
    SETERRQ(PETSC_COMM_SELF,MOFEM_DATA_INCONSISTENCY,"data inconsistency");
  }
  for(int ff = 0;ff<4;ff++) {
    data.dataOnEntities[MBTRI][ff].getHdivN(base).resize(G_DIM,3*NBFACETRI_HDIV_AINSWORTH_COLE(faces_order[ff]),false);
    data.dataOnEntities[MBTRI][ff].getDiffHdivN(base).resize(G_DIM,9*NBFACETRI_HDIV_AINSWORTH_COLE(faces_order[ff]),false);
    int col = 0,diff_col = 0;
    for(int oo = 0;oo<faces_order[ff];oo++) {
      for(int ee = 0;ee<3;ee++) {
        //values
        for(int dd = 3*NBFACETRI_EDGE_HDIV_AINSWORTH_COLE(oo);dd<3*NBFACETRI_EDGE_HDIV_AINSWORTH_COLE(oo+1);dd++,col++) {
          for(int gg = 0;gg<G_DIM;gg++) {
            data.dataOnEntities[MBTRI][ff].getHdivN(base)(gg,col) = N_face_edge(ff,ee)(gg,dd);
          }
        }
        //direvatives
        for(int dd = 9*NBFACETRI_EDGE_HDIV_AINSWORTH_COLE(oo);dd<9*NBFACETRI_EDGE_HDIV_AINSWORTH_COLE(oo+1);dd++,diff_col++) {
          for(int gg = 0;gg<G_DIM;gg++) {
            data.dataOnEntities[MBTRI][ff].getDiffHdivN(base)(gg,diff_col) = diffN_face_edge(ff,ee)(gg,dd);
          }
        }
      }
      //values
      for(int dd = 3*NBFACETRI_FACE_HDIV_AINSWORTH_COLE(oo);dd<3*NBFACETRI_FACE_HDIV_AINSWORTH_COLE(oo+1);dd++,col++) {
        for(int gg = 0;gg<G_DIM;gg++) {
          data.dataOnEntities[MBTRI][ff].getHdivN(base)(gg,col) = N_face_bubble[ff](gg,dd);
        }
      }
      //direvatives
      for(int dd = 9*NBFACETRI_FACE_HDIV_AINSWORTH_COLE(oo);dd<9*NBFACETRI_FACE_HDIV_AINSWORTH_COLE(oo+1);dd++,diff_col++) {
        for(int gg = 0;gg<G_DIM;gg++) {
          data.dataOnEntities[MBTRI][ff].getDiffHdivN(base)(gg,diff_col) = diffN_face_bubble[ff](gg,dd);
        }
      }
    }
  }

  //volume
  int col = 0,diff_col = 0;
  data.dataOnEntities[MBTET][0].getHdivN(base).resize(G_DIM,3*NBVOLUMETET_HDIV_AINSWORTH_COLE(volume_order),false);
  data.dataOnEntities[MBTET][0].getDiffHdivN(base).resize(G_DIM,9*NBVOLUMETET_HDIV_AINSWORTH_COLE(volume_order),false);
  for(int oo = 0;oo<volume_order;oo++) {
    for(int ee = 0;ee<6;ee++) {
      //values
      for(int dd = 3*NBVOLUMETET_EDGE_HDIV_AINSWORTH_COLE(oo);dd<3*NBVOLUMETET_EDGE_HDIV_AINSWORTH_COLE(oo+1);dd++,col++) {
        for(int gg = 0;gg<G_DIM;gg++) {
          data.dataOnEntities[MBTET][0].getHdivN(base)(gg,col) = N_volume_edge[ee](gg,dd);
        }
      }
      //direvatives
      for(int dd = 9*NBVOLUMETET_EDGE_HDIV_AINSWORTH_COLE(oo);dd<9*NBVOLUMETET_EDGE_HDIV_AINSWORTH_COLE(oo+1);dd++,diff_col++) {
        for(int gg = 0;gg<G_DIM;gg++) {
          data.dataOnEntities[MBTET][0].getDiffHdivN(base)(gg,diff_col) = diffN_volume_edge[ee](gg,dd);
        }
      }
    }
    for(int ff = 0;ff<4;ff++) {
      //values
      for(int dd = 3*NBVOLUMETET_FACE_HDIV_AINSWORTH_COLE(oo);dd<3*NBVOLUMETET_FACE_HDIV_AINSWORTH_COLE(oo+1);dd++,col++) {
        for(int gg = 0;gg<G_DIM;gg++) {
          data.dataOnEntities[MBTET][0].getHdivN(base)(gg,col) = N_volume_face[ff](gg,dd);
        }
      }
      //direvatives
      for(int dd = 9*NBVOLUMETET_FACE_HDIV_AINSWORTH_COLE(oo);dd<9*NBVOLUMETET_FACE_HDIV_AINSWORTH_COLE(oo+1);dd++,diff_col++) {
        for(int gg = 0;gg<G_DIM;gg++) {
          data.dataOnEntities[MBTET][0].getDiffHdivN(base)(gg,diff_col) = diffN_volume_face[ff](gg,dd);
        }
      }
    }
    //values
    for(int dd = 3*NBVOLUMETET_VOLUME_HDIV_AINSWORTH_COLE(oo);dd<3*NBVOLUMETET_VOLUME_HDIV_AINSWORTH_COLE(oo+1);dd++,col++) {
      for(int gg = 0;gg<G_DIM;gg++) {
        data.dataOnEntities[MBTET][0].getHdivN(base)(gg,col) = N_volume_bubble(gg,dd);
      }
    }
    //direvatives
    for(int dd = 9*NBVOLUMETET_VOLUME_HDIV_AINSWORTH_COLE(oo);dd<9*NBVOLUMETET_VOLUME_HDIV_AINSWORTH_COLE(oo+1);dd++,diff_col++) {
      for(int gg = 0;gg<G_DIM;gg++) {
        data.dataOnEntities[MBTET][0].getDiffHdivN(base)(gg,diff_col) = diffN_volume_bubble(gg,dd);
      }
    }

  }

  PetscFunctionReturn(0);
}

PetscErrorCode ForcesAndSurcesCore::shapeTRIFunctions_H1(
  DataForcesAndSurcesCore &data,
  const double *G_X,
  const double *G_Y,
  const int G_DIM,
  const FieldApproximationBase base,
  PetscErrorCode (*base_polynomials)(int p,double s,double *diff_s,double *L,double *diffL,const int dim)
) {
  PetscFunctionBegin;

  data.dataOnEntities[MBVERTEX][0].getDiffN(base).resize(3,2,false);
  ierr = ShapeDiffMBTRI(&*data.dataOnEntities[MBVERTEX][0].getDiffN(base).data().begin()); CHKERRQ(ierr);

  if(data.spacesOnEntities[MBEDGE].test(H1)) {
    //edges
    if(data.dataOnEntities[MBEDGE].size()!=3) {
      SETERRQ(PETSC_COMM_SELF,MOFEM_DATA_INCONSISTENCY,"data inconsistency");
    }
    int sense[3],order[3];
    double *H1edgeN[3],*diffH1edgeN[3];
    for(int ee = 0;ee<3;ee++) {
      if(data.dataOnEntities[MBEDGE][ee].getSense() == 0) {
        SETERRQ(PETSC_COMM_SELF,MOFEM_DATA_INCONSISTENCY,"data inconsistency");
      }
      sense[ee] = data.dataOnEntities[MBEDGE][ee].getSense();
      order[ee] = data.dataOnEntities[MBEDGE][ee].getDataOrder();
      int nb_dofs = NBEDGE_H1_AINSWORTH_COLE(data.dataOnEntities[MBEDGE][ee].getDataOrder());
      data.dataOnEntities[MBEDGE][ee].getN(base).resize(G_DIM,nb_dofs,false);
      data.dataOnEntities[MBEDGE][ee].getDiffN(base).resize(G_DIM,2*nb_dofs,false);
      H1edgeN[ee] = &*data.dataOnEntities[MBEDGE][ee].getN(base).data().begin();
      diffH1edgeN[ee] = &*data.dataOnEntities[MBEDGE][ee].getDiffN(base).data().begin();
    }
    ierr = H1_EdgeShapeFunctions_MBTRI(sense,order,
      &*data.dataOnEntities[MBVERTEX][0].getN(base).data().begin(),
      &*data.dataOnEntities[MBVERTEX][0].getDiffN(base).data().begin(),
      H1edgeN,diffH1edgeN,G_DIM,base_polynomials
    ); CHKERRQ(ierr);
  }

  if(data.spacesOnEntities[MBTRI].test(H1)) {
    //face
    if(data.dataOnEntities[MBTRI].size()!=1) {
      SETERRQ(PETSC_COMM_SELF,MOFEM_DATA_INCONSISTENCY,"data inconsistency");
    }
    int nb_dofs = NBFACETRI_H1_AINSWORTH_COLE(data.dataOnEntities[MBTRI][0].getDataOrder());
    data.dataOnEntities[MBTRI][0].getN(base).resize(G_DIM,nb_dofs,false);
    data.dataOnEntities[MBTRI][0].getDiffN(base).resize(G_DIM,2*nb_dofs,false);
    const int face_nodes[] = { 0,1,2 };
    ierr = H1_FaceShapeFunctions_MBTRI(
      face_nodes,
      data.dataOnEntities[MBTRI][0].getDataOrder(),
      &*data.dataOnEntities[MBVERTEX][0].getN(base).data().begin(),
      &*data.dataOnEntities[MBVERTEX][0].getDiffN(base).data().begin(),
      &*data.dataOnEntities[MBTRI][0].getN(base).data().begin(),
      &*data.dataOnEntities[MBTRI][0].getDiffN(base).data().begin(),
      G_DIM,base_polynomials
    ); CHKERRQ(ierr);
  }

  PetscFunctionReturn(0);
}

PetscErrorCode ForcesAndSurcesCore::shapeFlatPRISMFunctions_H1(
  DataForcesAndSurcesCore &data,
  const double *G_X,
  const double *G_Y,
  const int G_DIM,
  const FieldApproximationBase base,
  PetscErrorCode (*base_polynomials)(int p,double s,double *diff_s,double *L,double *diffL,const int dim)
) {
  PetscFunctionBegin;

  try {

    int num_nodes;
    const EntityHandle *conn_prism;
    rval = mField.get_moab().get_connectivity(fePtr->get_ent(),conn_prism,num_nodes,true); CHKERRQ_MOAB(rval);

    SideNumber_multiIndex& side_table = const_cast<SideNumber_multiIndex&>(fePtr->get_side_number_table());
    SideNumber_multiIndex::nth_index<1>::type::iterator siit3 = side_table.get<1>().find(boost::make_tuple(MBTRI,3));
    SideNumber_multiIndex::nth_index<1>::type::iterator siit4 = side_table.get<1>().find(boost::make_tuple(MBTRI,4));
    if(siit3==side_table.get<1>().end()) SETERRQ(PETSC_COMM_SELF,MOFEM_DATA_INCONSISTENCY,"data inconsistency");
    if(siit4==side_table.get<1>().end()) SETERRQ(PETSC_COMM_SELF,MOFEM_DATA_INCONSISTENCY,"data inconsistency");
    const EntityHandle *conn_face3;
    const EntityHandle *conn_face4;
    rval = mField.get_moab().get_connectivity(siit3->ent,conn_face3,num_nodes,true); CHKERRQ_MOAB(rval);
    rval = mField.get_moab().get_connectivity(siit4->ent,conn_face4,num_nodes,true); CHKERRQ_MOAB(rval);

    MatrixDouble N(G_DIM,3);
    ierr = ShapeMBTRI(&*N.data().begin(),G_X,G_Y,G_DIM); CHKERRQ(ierr);
    MatrixDouble diffN(3,2);
    ierr = ShapeDiffMBTRI(&*diffN.data().begin()); CHKERRQ(ierr);
    data.dataOnEntities[MBVERTEX][0].getN(base).resize(G_DIM,6,false);
    data.dataOnEntities[MBVERTEX][0].getDiffN(base).resize(G_DIM,12,false);

    int face_nodes[2][3];
    for(int nn = 0;nn<3;nn++) {
      face_nodes[0][nn] = distance(conn_prism,find(conn_prism,conn_prism+3,conn_face3[nn]));
      face_nodes[1][nn] = distance(conn_prism+3,find(conn_prism,conn_prism+6,conn_face4[nn]));
      // cerr << "fc " << 0 << " " << nn << " " << face_nodes[0][nn]  << endl;
      // cerr << "fc " << 1 << " " << nn << " " << face_nodes[1][nn]  << endl;
      for(int gg = 0;gg<G_DIM;gg++) {
        double val = N(gg,nn);
        double val_x = diffN(nn,0);
        double val_y = diffN(nn,1);
        data.dataOnEntities[MBVERTEX][0].getN(base)(gg,nn) = val;
        data.dataOnEntities[MBVERTEX][0].getDiffN(base)(gg,2*nn+0) = val_x;
        data.dataOnEntities[MBVERTEX][0].getDiffN(base)(gg,2*nn+1) = val_y;
        data.dataOnEntities[MBVERTEX][0].getN(base)(gg,3+nn) = val;
        data.dataOnEntities[MBVERTEX][0].getDiffN(base)(gg,6+2*nn+0) = val_x;
        data.dataOnEntities[MBVERTEX][0].getDiffN(base)(gg,6+2*nn+1) = val_y;
      }
    }

    //edges
    if(data.dataOnEntities[MBEDGE].size()!=9) {
      SETERRQ(PETSC_COMM_SELF,MOFEM_DATA_INCONSISTENCY,"data inconsistency");
    }
    int valid_edges[] = { 1,1,1, 0,0,0, 1,1,1 };
    int sense[9],order[9];
    double *H1edgeN[9],*diffH1edgeN[9];
    if((data.spacesOnEntities[MBEDGE]).test(H1)) {
      for(int ee = 0;ee<9;ee++) {
        if(!valid_edges[ee]) continue;
        if(data.dataOnEntities[MBEDGE][ee].getSense() == 0) {
          SETERRQ(PETSC_COMM_SELF,MOFEM_DATA_INCONSISTENCY,"data inconsistency");
        }
        sense[ee] = data.dataOnEntities[MBEDGE][ee].getSense();
        order[ee] = data.dataOnEntities[MBEDGE][ee].getDataOrder();
        int nb_dofs = NBEDGE_H1_AINSWORTH_COLE(data.dataOnEntities[MBEDGE][ee].getDataOrder());
        data.dataOnEntities[MBEDGE][ee].getN(base).resize(G_DIM,nb_dofs,false);
        data.dataOnEntities[MBEDGE][ee].getDiffN(base).resize(G_DIM,2*nb_dofs,false);
        H1edgeN[ee] = &*data.dataOnEntities[MBEDGE][ee].getN(base).data().begin();
        diffH1edgeN[ee] = &*data.dataOnEntities[MBEDGE][ee].getDiffN(base).data().begin();
      }
      //shape functions on face 3
      ierr = H1_EdgeShapeFunctions_MBTRI(
        &sense[0],
        &order[0],
        &*N.data().begin(),
        &*diffN.data().begin(),
        &H1edgeN[0],
        &diffH1edgeN[0],
        G_DIM,base_polynomials
      ); CHKERRQ(ierr);
      //shape functions on face 4
      ierr = H1_EdgeShapeFunctions_MBTRI(
        &sense[6],
        &order[6],
        &*N.data().begin(),
        &*diffN.data().begin(),
        &H1edgeN[6],
        &diffH1edgeN[6],
        G_DIM,base_polynomials
      ); CHKERRQ(ierr);
    }

    //face
    if(data.dataOnEntities[MBTRI].size()!=5) {
      SETERRQ(PETSC_COMM_SELF,MOFEM_DATA_INCONSISTENCY,"data inconsistency");
    }
    if((data.spacesOnEntities[MBTRI]).test(H1)) {
      for(int ff = 3;ff<=4;ff++) {
        int nb_dofs = NBFACETRI_H1_AINSWORTH_COLE(data.dataOnEntities[MBTRI][ff].getDataOrder());
        data.dataOnEntities[MBTRI][ff].getN(base).resize(G_DIM,nb_dofs,false);
        data.dataOnEntities[MBTRI][ff].getDiffN(base).resize(G_DIM,2*nb_dofs,false);
        ierr = H1_FaceShapeFunctions_MBTRI(
          face_nodes[ff-3],
          data.dataOnEntities[MBTRI][ff].getDataOrder(),
          &*N.data().begin(),
          &*diffN.data().begin(),
          &*data.dataOnEntities[MBTRI][ff].getN(base).data().begin(),
          &*data.dataOnEntities[MBTRI][ff].getDiffN(base).data().begin(),
          G_DIM,base_polynomials
        ); CHKERRQ(ierr);
      }
    }

    } catch (MoFEMException const &e) {
      SETERRQ(PETSC_COMM_SELF,e.errorCode,e.errorMessage);
    } catch (exception& ex) {
      ostringstream ss;
      ss << "thorw in method: " << ex.what() << " at line " << __LINE__ << " in file " << __FILE__;
      SETERRQ(PETSC_COMM_SELF,MOFEM_STD_EXCEPTION_THROW,ss.str().c_str());
    }

  PetscFunctionReturn(0);
}

PetscErrorCode ForcesAndSurcesCore::shapeTRIFunctions_Hdiv(
  DataForcesAndSurcesCore &data,
  const double *G_X,
  const double *G_Y,
  const int G_DIM,
  const FieldApproximationBase base,
  PetscErrorCode (*base_polynomials)(int p,double s,double *diff_s,double *L,double *diffL,const int dim)
) {
  PetscFunctionBegin;

  data.dataOnEntities[MBVERTEX][0].getN(base).resize(G_DIM,3,false);
  ierr = ShapeMBTRI(&*data.dataOnEntities[MBVERTEX][0].getN(base).data().begin(),G_X,G_Y,G_DIM); CHKERRQ(ierr);

  double *PHI_f_e[3];
  double *PHI_f;

  N_face_edge.resize(1,3,false);
  N_face_bubble.resize(1,false);
  int face_order = data.dataOnEntities[MBTRI][0].getDataOrder();
  //three edges on face
  for(int ee = 0;ee<3;ee++) {
    N_face_edge(0,ee).resize(G_DIM,3*NBFACETRI_EDGE_HDIV_AINSWORTH_COLE(face_order),false);
    PHI_f_e[ee] = &((N_face_edge(0,ee))(0,0));
  }
  N_face_bubble[0].resize(G_DIM,3*NBFACETRI_FACE_HDIV_AINSWORTH_COLE(face_order),false);
  PHI_f = &*(N_face_bubble[0].data().begin());

  int face_nodes[3] = { 0,1,2 };
  ierr = Hdiv_EdgeFaceShapeFunctions_MBTET_ON_FACE(
    face_nodes,face_order,
    &data.dataOnEntities[MBVERTEX][0].getN(base)(0,0),
    NULL,
    PHI_f_e,NULL,G_DIM,3,
    base_polynomials
  ); CHKERRQ(ierr);
  ierr = Hdiv_FaceBubbleShapeFunctions_MBTET_ON_FACE(
    face_nodes,face_order,
    &data.dataOnEntities[MBVERTEX][0].getN(base)(0,0),NULL,
    PHI_f,NULL,G_DIM,3,
    base_polynomials
  ); CHKERRQ(ierr);

  // set shape functions into data structure

  if(data.dataOnEntities[MBTRI].size()!=1) {
    SETERRQ(PETSC_COMM_SELF,MOFEM_DATA_INCONSISTENCY,"data inconsistency");
  }
  data.dataOnEntities[MBTRI][0].getHdivN(base).resize(G_DIM,3*NBFACETRI_HDIV_AINSWORTH_COLE(face_order),false);
  int col = 0;
  for(int oo = 0;oo<face_order;oo++) {
    for(int ee = 0;ee<3;ee++) {
      for(int dd = 3*NBFACETRI_EDGE_HDIV_AINSWORTH_COLE(oo);dd<3*NBFACETRI_EDGE_HDIV_AINSWORTH_COLE(oo+1);dd++,col++) {
        for(int gg = 0;gg<G_DIM;gg++) {
          data.dataOnEntities[MBTRI][0].getHdivN(base)(gg,col) = N_face_edge(0,ee)(gg,dd);
        }
      }
    }
    for(int dd = 3*NBFACETRI_FACE_HDIV_AINSWORTH_COLE(oo);dd<3*NBFACETRI_FACE_HDIV_AINSWORTH_COLE(oo+1);dd++,col++) {
      for(int gg = 0;gg<G_DIM;gg++) {
        data.dataOnEntities[MBTRI][0].getHdivN(base)(gg,col) = N_face_bubble[0](gg,dd);
      }
    }
  }

  PetscFunctionReturn(0);
}

PetscErrorCode ForcesAndSurcesCore::shapeEDGEFunctions_H1(
  DataForcesAndSurcesCore &data,
  int side_number,
  const double *G_X,
  const int G_DIM,
  const FieldApproximationBase base,
  PetscErrorCode (*base_polynomials)(int p,double s,double *diff_s,double *L,double *diffL,const int dim)
) {
  PetscFunctionBegin;

  // data.dataOnEntities[MBVERTEX][0].getN(base).resize(G_DIM,2,false);
  // ierr = ShapeMBEDGE(
  //   &*data.dataOnEntities[MBVERTEX][0].getN(base).data().begin(),
  //   G_X,
  //   G_DIM
  // ); CHKERRQ(ierr);
  data.dataOnEntities[MBVERTEX][0].getDiffN(base).resize(2,1,false);
  ierr = ShapeDiffMBEDGE(
    &*data.dataOnEntities[MBVERTEX][0].getDiffN(base).data().begin()
  ); CHKERRQ(ierr);

  //cerr << data.dataOnEntities[MBVERTEX][0].getN(base) << endl;
  //cerr << data.dataOnEntities[MBVERTEX][0].getDiffN(base) << endl;

  int sense = data.dataOnEntities[MBEDGE][side_number].getSense();
  int order = data.dataOnEntities[MBEDGE][side_number].getDataOrder();
  data.dataOnEntities[MBEDGE][side_number].getN(base).resize(G_DIM,NBEDGE_H1_AINSWORTH_COLE(order),false);
  data.dataOnEntities[MBEDGE][side_number].getDiffN(base).resize(G_DIM,NBEDGE_H1_AINSWORTH_COLE(order),false);
  if(data.dataOnEntities[MBEDGE][side_number].getDataOrder()>1) {
    double diff_s = 2.; // s = s(xi), ds/dxi = 2., because change of basis
    for(int gg = 0;gg<G_DIM;gg++) {
      double s = 2*G_X[gg]-1; // makes form -1..1
      if(!sense) {
        s *= -1;
        diff_s *= -1;
      }

      // calculate Legendre polynomials at integration points
      ierr = base_polynomials(
        NBEDGE_H1_AINSWORTH_COLE(order)-1,
        s,
        &diff_s,
        &data.dataOnEntities[MBEDGE][side_number].getN(base)(gg,0),
        &data.dataOnEntities[MBEDGE][side_number].getDiffN(base)(gg,0),1
      ); CHKERRQ(ierr);

      for(unsigned int pp = 0;pp<data.dataOnEntities[MBEDGE][side_number].getN(base).size2();pp++) {

        double L = data.dataOnEntities[MBEDGE][side_number].getN(base)(gg,pp);
        double diffL = data.dataOnEntities[MBEDGE][side_number].getDiffN(base)(gg,pp);

        // Calculate edge shape functions N0*N1*L(p), where N0 and N1 are nodal shape functions
        data.dataOnEntities[MBEDGE][side_number].getN(base)(gg,pp) =
          data.dataOnEntities[MBVERTEX][0].getN(base)(gg,0)*
          data.dataOnEntities[MBVERTEX][0].getN(base)(gg,1)*L;

        // Calculate derivative edge shape functions
        // dN/dksi = dN0/dxi*N1*L + N0*dN1/ksi*L + N0*N1*dL/dxi
        data.dataOnEntities[MBEDGE][side_number].getDiffN(base)(gg,pp) =
          ((+1.)*data.dataOnEntities[MBVERTEX][0].getN(base)(gg,1)
          +data.dataOnEntities[MBVERTEX][0].getN(base)(gg,0)*(-1.))*L
          +data.dataOnEntities[MBVERTEX][0].getN(base)(gg,0)*
          data.dataOnEntities[MBVERTEX][0].getN(base)(gg,1)*diffL;
      }
    }
  }

  //cerr << data.dataOnEntities[MBEDGE][0].getN(base) << endl;
  //cerr << data.dataOnEntities[MBEDGE][0].getDiffN(base) << endl;

  PetscFunctionReturn(0);
}

PetscErrorCode ForcesAndSurcesCore::shapeFlatPRISMFunctions_Hdiv(
  DataForcesAndSurcesCore &data,
  const double *G_X,
  const double *G_Y,
  const int G_DIM,
  const FieldApproximationBase base,
  PetscErrorCode (*base_polynomials)(int p,double s,double *diff_s,double *L,double *diffL,const int dim)
) {
  PetscFunctionBegin;
  SETERRQ(PETSC_COMM_SELF,MOFEM_NOT_IMPLEMENTED,"not yet implemented, i.e. flat prism with Hdiv space");
  PetscFunctionReturn(0);
}

// **** Data Operator ****

static PetscErrorCode get_porblem_row_indices(
  const ForcesAndSurcesCore *fe_ptr,const EntityType type,const int side,const string field_name,VectorInt& indices) {
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
  const ForcesAndSurcesCore *fe_ptr,const EntityType type,const int side,const string field_name,VectorInt& indices) {
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
  const string field_name,const EntityType type,const int side,VectorInt& indices) {
  PetscFunctionBegin;

  PetscErrorCode ierr;

  if(ptrFE == NULL) {
    SETERRQ(PETSC_COMM_SELF,MOFEM_DATA_INCONSISTENCY,"data inconsistency");
  }
  ierr = get_porblem_row_indices(ptrFE,type,side,field_name,indices); CHKERRQ(ierr);

  PetscFunctionReturn(0);
}

PetscErrorCode ForcesAndSurcesCore::UserDataOperator::getPorblemColIndices(
  const string field_name,const EntityType type,const int side,VectorInt& indices) {
  PetscFunctionBegin;

  PetscErrorCode ierr;

  if(ptrFE == NULL) {
    SETERRQ(PETSC_COMM_SELF,MOFEM_DATA_INCONSISTENCY,"data inconsistency");
  }
  ierr = get_porblem_col_indices(ptrFE,type,side,field_name,indices); CHKERRQ(ierr);

  PetscFunctionReturn(0);
}

}
