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

#include <petscsys.h>
#include <petscvec.h>
#include <petscmat.h>
#include <petscsnes.h>
#include <petscts.h>

#include <definitions.h>
#include <h1_hdiv_hcurl_l2.h>
#include <fem_tools.h>

#include <Common.hpp>

#include <LoopMethods.hpp>
#include <FieldInterface.hpp>

#define BOOST_UBLAS_SHALLOW_ARRAY_ADAPTOR

#include <boost/numeric/ublas/storage.hpp>
#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/matrix_proxy.hpp>
#include <boost/numeric/ublas/io.hpp>

#include <boost/ptr_container/ptr_vector.hpp>
#include <boost/ptr_container/ptr_map.hpp>

#include <ForcesAndSurcesCore.hpp>

#ifdef __cplusplus
extern "C" {
#endif
  #include <cblas.h>
  #include <lapack_wrap.h>
  #include <gm_rule.h>
#ifdef __cplusplus
}
#endif

static ErrorCode rval;

namespace MoFEM {

// ** Sense **

PetscErrorCode ForcesAndSurcesCore::getSense(EntityType type,boost::ptr_vector<DataForcesAndSurcesCore::EntData> &data) {
  PetscFunctionBegin;
  try {
    SideNumber_multiIndex& side_table = const_cast<SideNumber_multiIndex&>(fePtr->get_side_number_table());
    if(data.size() < side_table.get<2>().count(type)) {
      // prims has 9 edges, some of edges for "flat" prism are not active
      SETERRQ2(PETSC_COMM_SELF,MOFEM_DATA_INCONSISTENCT,"data inconsistency %u < %u",data.size(),side_table.get<2>().count(type));
    }
    SideNumber_multiIndex::nth_index<2>::type::iterator siit = side_table.get<2>().lower_bound(type);
    SideNumber_multiIndex::nth_index<2>::type::iterator hi_siit = side_table.get<2>().upper_bound(type);
    for(;siit!=hi_siit;siit++) {
      data[siit->side_number].getSense() = siit->sense;
      if(siit->brother_side_number!=-1) {
        if(data.size() < (unsigned)siit->brother_side_number) {
          SETERRQ(PETSC_COMM_SELF,MOFEM_DATA_INCONSISTENCT,"data inconsistency");
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

// ** Order **

PetscErrorCode ForcesAndSurcesCore::getOrder(const EntityType type,const FieldSpace space,boost::ptr_vector<DataForcesAndSurcesCore::EntData> &data) {
  PetscFunctionBegin;
  try {
    SideNumber_multiIndex& side_table = const_cast<SideNumber_multiIndex&>(fePtr->get_side_number_table());
    if(data.size() < side_table.get<2>().count(type)) { // prims has 9 edges, some of edges for "flat" prism are not active
    SETERRQ2(PETSC_COMM_SELF,MOFEM_DATA_INCONSISTENCT,
      "data inconsistency %d < %d",
      data.size(),side_table.get<2>().count(type));
    }
    for(unsigned int side = 0;side<data.size();side++) {
      data[side].getOrder() = 0;
    }
    FEDofMoFEMEntity_multiIndex::index<Composite_EntType_and_Space_mi_tag>::type &data_dofs =
    const_cast<FEDofMoFEMEntity_multiIndex::index<Composite_EntType_and_Space_mi_tag>::type&>(fePtr->get_data_dofs().get<Composite_EntType_and_Space_mi_tag>());
    FEDofMoFEMEntity_multiIndex::index<Composite_EntType_and_Space_mi_tag>::type::iterator dit,hi_dit;
    dit = data_dofs.lower_bound(boost::make_tuple(type,space));
    hi_dit = data_dofs.upper_bound(boost::make_tuple(type,space));
    for(;dit!=hi_dit;dit++) {
      ApproximationOrder ent_order = dit->get_max_order();
      int side_number = dit->side_number_ptr->side_number;
      if(side_number < 0) {
        SETERRQ(PETSC_COMM_SELF,MOFEM_DATA_INCONSISTENCT,"data inconsistency");
      }
      data[side_number].getOrder() = data[side_number].getOrder() > ent_order ? data[side_number].getOrder() : ent_order;
      if(dit->side_number_ptr->brother_side_number!=-1) {
        if(data.size() < (unsigned int)dit->side_number_ptr->brother_side_number) {
          SETERRQ(PETSC_COMM_SELF,MOFEM_DATA_INCONSISTENCT,"data inconsistency");
        }
        data[dit->side_number_ptr->brother_side_number].getOrder() = data[side_number].getOrder();
      }
    }
  } catch (exception& ex) {
    ostringstream ss;
    ss << "thorw in method: " << ex.what() << " at line " << __LINE__ << " in file " << __FILE__;
    SETERRQ(PETSC_COMM_SELF,MOFEM_STD_EXCEPTION_THROW,ss.str().c_str());
  }
  PetscFunctionReturn(0);
}


PetscErrorCode ForcesAndSurcesCore::getEdgesOrder(DataForcesAndSurcesCore &data,const FieldSpace space) {
  PetscFunctionBegin;
  ierr = getOrder(MBEDGE,space,data.dataOnEntities[MBEDGE]); CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

PetscErrorCode ForcesAndSurcesCore::getTrisOrder(DataForcesAndSurcesCore &data,const FieldSpace space) {
  PetscFunctionBegin;
  ierr = getOrder(MBTRI,space,data.dataOnEntities[MBTRI]); CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

PetscErrorCode ForcesAndSurcesCore::getTetsOrder(DataForcesAndSurcesCore &data,const FieldSpace space) {
  PetscFunctionBegin;
  ierr = getOrder(MBTET,space,data.dataOnEntities[MBTET]); CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

PetscErrorCode ForcesAndSurcesCore::getOrder(const string &field_name,const EntityType type,boost::ptr_vector<DataForcesAndSurcesCore::EntData> &data) {
  PetscFunctionBegin;
  SideNumber_multiIndex& side_table = const_cast<SideNumber_multiIndex&>(fePtr->get_side_number_table());
  if(data.size() < side_table.get<2>().count(type)) { // prims has 9 edges, some of edges for "flat" prism are not active
  SETERRQ2(PETSC_COMM_SELF,MOFEM_DATA_INCONSISTENCT,
    "data inconsistency %d < %d",
    data.size(),side_table.get<2>().count(type));
  }
  for(unsigned int side = 0;side<data.size();side++) {
    data[side].getOrder() = 0;
  }
  FEDofMoFEMEntity_multiIndex::index<Composite_Name_And_Type_mi_tag>::type &data_dofs =
  const_cast<FEDofMoFEMEntity_multiIndex::index<Composite_Name_And_Type_mi_tag>::type&>(fePtr->get_data_dofs().get<Composite_Name_And_Type_mi_tag>());
  FEDofMoFEMEntity_multiIndex::index<Composite_Name_And_Type_mi_tag>::type::iterator dit,hi_dit;
  dit = data_dofs.lower_bound(boost::make_tuple(field_name,type));
  hi_dit = data_dofs.upper_bound(boost::make_tuple(field_name,type));
  for(;dit!=hi_dit;dit++) {
    ApproximationOrder ent_order = dit->get_max_order();
    int side_number = dit->side_number_ptr->side_number;
    if(side_number < 0) {
      SETERRQ(PETSC_COMM_SELF,MOFEM_DATA_INCONSISTENCT,"data inconsistency");
    }
    data[side_number].getOrder() = data[side_number].getOrder() > ent_order ? data[side_number].getOrder() : ent_order;
    if(dit->side_number_ptr->brother_side_number!=-1) {
      if(data.size() < (unsigned int)dit->side_number_ptr->brother_side_number) {
        SETERRQ(PETSC_COMM_SELF,MOFEM_DATA_INCONSISTENCT,"data inconsistency");
      }
      data[dit->side_number_ptr->brother_side_number].getOrder() = data[side_number].getOrder();
    }
  }
  PetscFunctionReturn(0);
}

PetscErrorCode ForcesAndSurcesCore::getEdgesOrder(DataForcesAndSurcesCore &data,const string &field_name) {
  PetscFunctionBegin;
  ierr = getOrder(field_name,MBEDGE,data.dataOnEntities[MBEDGE]); CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

PetscErrorCode ForcesAndSurcesCore::getTrisOrder(DataForcesAndSurcesCore &data,const string &field_name) {
  PetscFunctionBegin;
  ierr = getOrder(field_name,MBTRI,data.dataOnEntities[MBTRI]); CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

PetscErrorCode ForcesAndSurcesCore::getTetsOrder(DataForcesAndSurcesCore &data,const string &field_name) {
  PetscFunctionBegin;
  ierr = getOrder(field_name,MBTET,data.dataOnEntities[MBTET]); CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

// ** Indices **

PetscErrorCode ForcesAndSurcesCore::getNodesIndices(
  const string &field_name,FENumeredDofMoFEMEntity_multiIndex &dofs,ublas::vector<int> &nodes_indices
) {
  PetscFunctionBegin;
  FENumeredDofMoFEMEntity_multiIndex::index<Composite_Name_And_Type_mi_tag>::type::iterator dit,hi_dit;
  dit = dofs.get<Composite_Name_And_Type_mi_tag>().lower_bound(boost::make_tuple(field_name,MBVERTEX));
  hi_dit = dofs.get<Composite_Name_And_Type_mi_tag>().upper_bound(boost::make_tuple(field_name,MBVERTEX));
  nodes_indices.resize(distance(dit,hi_dit));
  for(;dit!=hi_dit;dit++) {
    int idx = dit->get_petsc_gloabl_dof_idx();
    int side_number = dit->side_number_ptr->side_number;
    nodes_indices[side_number*dit->get_max_rank()+dit->get_dof_rank()] = idx;
    int  brother_side_number = dit->side_number_ptr->brother_side_number;
    if(brother_side_number!=-1) {
      if(nodes_indices.size()<(unsigned int)(brother_side_number*dit->get_max_rank()+dit->get_max_rank())) {
        nodes_indices.resize(brother_side_number*dit->get_max_rank()+dit->get_max_rank());
      }
      nodes_indices[brother_side_number*dit->get_max_rank()+dit->get_dof_rank()] = idx;
    }
  }
  PetscFunctionReturn(0);
}

PetscErrorCode ForcesAndSurcesCore::getRowNodesIndices(DataForcesAndSurcesCore &data,const string &field_name) {
  PetscFunctionBegin;
  ierr = getNodesIndices(
    field_name,const_cast<FENumeredDofMoFEMEntity_multiIndex&>(fePtr->get_rows_dofs()),data.dataOnEntities[MBVERTEX][0].getIndices()
  ); CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

PetscErrorCode ForcesAndSurcesCore::getColNodesIndices(DataForcesAndSurcesCore &data,const string &field_name) {
  PetscFunctionBegin;
  ierr = getNodesIndices(field_name,const_cast<FENumeredDofMoFEMEntity_multiIndex&>(fePtr->get_cols_dofs()),data.dataOnEntities[MBVERTEX][0].getIndices()); CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

PetscErrorCode ForcesAndSurcesCore::getTypeIndices(
  const string &field_name,FENumeredDofMoFEMEntity_multiIndex &dofs,EntityType type,int side_number,ublas::vector<int> &indices
) {
  PetscFunctionBegin;
  FENumeredDofMoFEMEntity_multiIndex::index<Composite_Name_Type_And_Side_Number_mi_tag>::type::iterator dit,hi_dit;
  dit = dofs.get<Composite_Name_Type_And_Side_Number_mi_tag>().lower_bound(boost::make_tuple(field_name,type,side_number));
  hi_dit = dofs.get<Composite_Name_Type_And_Side_Number_mi_tag>().upper_bound(boost::make_tuple(field_name,type,side_number));
  indices.resize(0);
  for(;dit!=hi_dit;dit++) {
    int idx = dit->get_petsc_gloabl_dof_idx();
    indices.resize(dit->get_nb_dofs_on_ent());
    indices[dit->get_EntDofIdx()] = idx;
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
    ierr = getTypeIndices(field_name,dofs,type,siit->side_number,data[siit->side_number].getIndices()); CHKERRQ(ierr);
    if(siit->brother_side_number!=-1) {
      ierr = getTypeIndices(field_name,dofs,type,siit->side_number,data[siit->brother_side_number].getIndices()); CHKERRQ(ierr);
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
    SETERRQ(PETSC_COMM_SELF,MOFEM_DATA_INCONSISTENCT,"data inconsistency");
  }
  ierr = getTypeIndices(
    field_name,const_cast<FENumeredDofMoFEMEntity_multiIndex&>(fePtr->get_rows_dofs()),MBTET,0,data.dataOnEntities[MBTET][0].getIndices()
  ); CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

PetscErrorCode ForcesAndSurcesCore::getTetsColIndices(DataForcesAndSurcesCore &data,const string &field_name) {
  PetscFunctionBegin;
  if(data.dataOnEntities[MBTET].size() == 0) {
    SETERRQ(PETSC_COMM_SELF,MOFEM_DATA_INCONSISTENCT,"data inconsistency");
  }
  ierr = getTypeIndices(
    field_name,const_cast<FENumeredDofMoFEMEntity_multiIndex&>(fePtr->get_cols_dofs()),MBTET,0,data.dataOnEntities[MBTET][0].getIndices()
  ); CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

PetscErrorCode ForcesAndSurcesCore::getNoFieldIndices(
  const string &field_name,FENumeredDofMoFEMEntity_multiIndex &dofs,ublas::vector<int> &nodes_indices
) {
  PetscFunctionBegin;
  FENumeredDofMoFEMEntity_multiIndex::index<FieldName_mi_tag>::type::iterator dit,hi_dit;
  dit = dofs.get<FieldName_mi_tag>().lower_bound(field_name);
  hi_dit = dofs.get<FieldName_mi_tag>().upper_bound(field_name);
  nodes_indices.resize(distance(dit,hi_dit));
  for(;dit!=hi_dit;dit++) {
    int idx = dit->get_petsc_gloabl_dof_idx();
    nodes_indices[dit->get_dof_rank()] = idx;
  }
  PetscFunctionReturn(0);
}


PetscErrorCode ForcesAndSurcesCore::getNoFieldRowIndices(DataForcesAndSurcesCore &data,const string &field_name) {
  PetscFunctionBegin;
  //EntityType fe_type = fePtr->get_ent_type();
  if(data.dataOnEntities[MBENTITYSET].size() == 0) {
    SETERRQ(PETSC_COMM_SELF,MOFEM_DATA_INCONSISTENCT,"data inconsistency");
  }
  ierr = getNoFieldIndices(
    field_name,const_cast<FENumeredDofMoFEMEntity_multiIndex&>(fePtr->get_rows_dofs()),data.dataOnEntities[MBENTITYSET][0].getIndices()
  ); CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

PetscErrorCode ForcesAndSurcesCore::getNoFieldColIndices(DataForcesAndSurcesCore &data,const string &field_name) {
  PetscFunctionBegin;
  if(data.dataOnEntities[MBENTITYSET].size() == 0) {
    SETERRQ(PETSC_COMM_SELF,MOFEM_DATA_INCONSISTENCT,"data inconsistency");
  }
  ierr = getNoFieldIndices(
    field_name,const_cast<FENumeredDofMoFEMEntity_multiIndex&>(fePtr->get_cols_dofs()),data.dataOnEntities[MBENTITYSET][0].getIndices()
  ); CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

// ** Indices from problem **

PetscErrorCode ForcesAndSurcesCore::getProblemNodesIndices(const string &field_name,
  const NumeredDofMoFEMEntity_multiIndex &dofs,
  ublas::vector<int> &nodes_indices
) const {
  PetscFunctionBegin;

  nodes_indices.resize(0);

  SideNumber_multiIndex& side_table = const_cast<SideNumber_multiIndex&>(fePtr->get_side_number_table());
  SideNumber_multiIndex::nth_index<2>::type::iterator siit = side_table.get<2>().lower_bound(MBVERTEX);
  SideNumber_multiIndex::nth_index<2>::type::iterator hi_siit = side_table.get<2>().upper_bound(MBVERTEX);

  int nn = 0;
  for(;siit!=hi_siit;siit++,nn++) {

    const EntityHandle ent = siit->ent;
    NumeredDofMoFEMEntity_multiIndex::index<Composite_Name_And_Ent_And_EndDofIdx_mi_tag>::type::iterator dit,hi_dit;
    dit = dofs.get<Composite_Name_And_Ent_And_EndDofIdx_mi_tag>().lower_bound(boost::make_tuple(field_name,ent,0));
    hi_dit = dofs.get<Composite_Name_And_Ent_And_EndDofIdx_mi_tag>().upper_bound(boost::make_tuple(field_name,ent,10000));  /// very large number

    if(dit!=hi_dit) {

      if(!nn) {
        nodes_indices.resize(dit->get_max_rank()*distance(siit,hi_siit));
      }
      for(;dit!=hi_dit;dit++) {
        nodes_indices[siit->side_number*dit->get_max_rank()+dit->get_dof_rank()] = dit->get_petsc_gloabl_dof_idx();
      }

    }

  }

  PetscFunctionReturn(0);

}

PetscErrorCode ForcesAndSurcesCore::getProblemTypeIndices(
  const string &field_name,const NumeredDofMoFEMEntity_multiIndex &dofs,
  EntityType type,int side_number,ublas::vector<int> &indices
) const {
  PetscFunctionBegin;

  indices.resize(0);

  SideNumber_multiIndex& side_table = const_cast<SideNumber_multiIndex&>(fePtr->get_side_number_table());
  SideNumber_multiIndex::nth_index<1>::type::iterator siit = side_table.get<1>().lower_bound(boost::make_tuple(type,side_number));
  SideNumber_multiIndex::nth_index<1>::type::iterator hi_siit = side_table.get<1>().upper_bound(boost::make_tuple(type,side_number));

  for(;siit!=hi_siit;siit++) {

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

PetscErrorCode ForcesAndSurcesCore::getProblemNodesRowIndices(const string &field_name,ublas::vector<int> &nodes_indices) const {
  PetscFunctionBegin;
  PetscErrorCode ierr;
  ierr = getProblemNodesIndices(field_name,problemPtr->numered_dofs_rows,nodes_indices); CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

PetscErrorCode ForcesAndSurcesCore::getProblemTypeRowIndices(const string &field_name,EntityType type,int side_number,ublas::vector<int> &indices) const {
  PetscFunctionBegin;
  PetscErrorCode ierr;
  ierr = getProblemTypeIndices(field_name,problemPtr->numered_dofs_rows,type,side_number,indices); CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

PetscErrorCode ForcesAndSurcesCore::getProblemNodesColIndices(const string &field_name,ublas::vector<int> &nodes_indices) const {
  PetscFunctionBegin;
  PetscErrorCode ierr;
  ierr = getProblemNodesIndices(field_name,problemPtr->numered_dofs_cols,nodes_indices); CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

PetscErrorCode ForcesAndSurcesCore::getProblemTypeColIndices(const string &field_name,EntityType type,int side_number,ublas::vector<int> &indices) const {
  PetscFunctionBegin;
  PetscErrorCode ierr;
  ierr = getProblemTypeIndices(field_name,problemPtr->numered_dofs_cols,type,side_number,indices); CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

// ** Data **

PetscErrorCode ForcesAndSurcesCore::getNodesFieldData(
  const string &field_name,FEDofMoFEMEntity_multiIndex &dofs,ublas::vector<double> &nodes_field_data
) {
  PetscFunctionBegin;
  try {
    FEDofMoFEMEntity_multiIndex::index<Composite_Name_And_Type_mi_tag>::type::iterator dit,hi_dit;
    dit = dofs.get<Composite_Name_And_Type_mi_tag>().lower_bound(boost::make_tuple(field_name,MBVERTEX));
    hi_dit = dofs.get<Composite_Name_And_Type_mi_tag>().upper_bound(boost::make_tuple(field_name,MBVERTEX));
    nodes_field_data.resize(distance(dit,hi_dit));
    for(;dit!=hi_dit;dit++) {
      FieldData val = dit->get_FieldData();
      int side_number = dit->side_number_ptr->side_number;
      if(side_number == -1) {
        SETERRQ(PETSC_COMM_SELF,MOFEM_DATA_INCONSISTENCT,"data inconsistency");
      }
      nodes_field_data[side_number*dit->get_max_rank()+dit->get_dof_rank()] = val;
      int  brother_side_number = dit->side_number_ptr->brother_side_number;
      if(brother_side_number!=-1) {
        if(nodes_field_data.size()<(unsigned int)(brother_side_number*dit->get_max_rank()+dit->get_max_rank())) {
          nodes_field_data.resize(brother_side_number*dit->get_max_rank()+dit->get_max_rank());
        }
        nodes_field_data[brother_side_number*dit->get_max_rank()+dit->get_dof_rank()] = val;
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
    field_name,const_cast<FEDofMoFEMEntity_multiIndex&>(fePtr->get_data_dofs()),data.dataOnEntities[MBVERTEX][0].getFieldData()
  ); CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

PetscErrorCode ForcesAndSurcesCore::getTypeFieldData(
  const string &field_name,FEDofMoFEMEntity_multiIndex &dofs,EntityType type,int side_number,ublas::vector<double> &ent_field_data
) {
  PetscFunctionBegin;
  FEDofMoFEMEntity_multiIndex::index<Composite_Name_Type_And_Side_Number_mi_tag>::type::iterator dit,hi_dit;
  dit = dofs.get<Composite_Name_Type_And_Side_Number_mi_tag>().lower_bound(boost::make_tuple(field_name,type,side_number));
  hi_dit = dofs.get<Composite_Name_Type_And_Side_Number_mi_tag>().upper_bound(boost::make_tuple(field_name,type,side_number));
  ent_field_data.resize(0);
  for(;dit!=hi_dit;dit++) {
    FieldData val = dit->get_FieldData();
    ent_field_data.resize(dit->get_nb_dofs_on_ent());
    ent_field_data[dit->get_EntDofIdx()] = val;
  }
  PetscFunctionReturn(0);
}

PetscErrorCode ForcesAndSurcesCore::getTypeFieldData(
  const string &field_name,FEDofMoFEMEntity_multiIndex &dofs,EntityType type,
  boost::ptr_vector<DataForcesAndSurcesCore::EntData> &data
) {
  PetscFunctionBegin;
  SideNumber_multiIndex& side_table = const_cast<SideNumber_multiIndex&>(fePtr->get_side_number_table());
  if(data.size() < side_table.get<2>().count(type)) {
    SETERRQ(PETSC_COMM_SELF,MOFEM_DATA_INCONSISTENCT,"data inconsistency");
  }
  SideNumber_multiIndex::nth_index<2>::type::iterator siit = side_table.get<2>().lower_bound(type);
  SideNumber_multiIndex::nth_index<2>::type::iterator hi_siit = side_table.get<2>().upper_bound(type);
  for(;siit!=hi_siit;siit++) {
    ierr = getTypeFieldData(field_name,dofs,type,siit->side_number,data[siit->side_number].getFieldData()); CHKERRQ(ierr);
    if(siit->brother_side_number!=-1) {
      if(data.size() < (unsigned int)siit->brother_side_number) {
        SETERRQ(PETSC_COMM_SELF,MOFEM_DATA_INCONSISTENCT,"data inconsistency");
      }
      ierr = getTypeFieldData(field_name,dofs,type,siit->side_number,data[siit->brother_side_number].getFieldData()); CHKERRQ(ierr);
    }
  }
  PetscFunctionReturn(0);
}

PetscErrorCode ForcesAndSurcesCore::getNoFieldFieldData(
  const string &field_name,FEDofMoFEMEntity_multiIndex &dofs,ublas::vector<double> &ent_field_data
) {
  PetscFunctionBegin;
  FEDofMoFEMEntity_multiIndex::index<FieldName_mi_tag>::type::iterator dit,hi_dit;
  dit = dofs.get<FieldName_mi_tag>().lower_bound(field_name);
  hi_dit = dofs.get<FieldName_mi_tag>().upper_bound(field_name);
  ent_field_data.resize(distance(dit,hi_dit));
  for(;dit!=hi_dit;dit++) {
    ent_field_data[dit->get_dof_rank()] = dit->get_FieldData();
  }
  PetscFunctionReturn(0);
}

PetscErrorCode ForcesAndSurcesCore::getNoFieldFieldData(
  DataForcesAndSurcesCore &data,const string &field_name
) {
  PetscFunctionBegin;
  if(data.dataOnEntities[MBENTITYSET].size() == 0) {
    SETERRQ(PETSC_COMM_SELF,MOFEM_DATA_INCONSISTENCT,"data inconsistency");
  }
  ierr = getNoFieldFieldData(
    field_name,const_cast<FEDofMoFEMEntity_multiIndex&>(fePtr->get_data_dofs()),data.dataOnEntities[MBENTITYSET][0].getFieldData()
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
    DataForcesAndSurcesCore &data,const string &field_name) {
    PetscFunctionBegin;
    ierr = getTypeFieldData(field_name,const_cast<FEDofMoFEMEntity_multiIndex&>(fePtr->get_data_dofs()),MBTRI,data.dataOnEntities[MBTRI]); CHKERRQ(ierr);
    PetscFunctionReturn(0);
  }

PetscErrorCode ForcesAndSurcesCore::getTetsFieldData(DataForcesAndSurcesCore &data,const string &field_name) {
  PetscFunctionBegin;
  if(data.dataOnEntities[MBTET].size() == 0) {
    SETERRQ(PETSC_COMM_SELF,MOFEM_DATA_INCONSISTENCT,"data inconsistency");
  }
  ierr = getTypeFieldData(field_name,const_cast<FEDofMoFEMEntity_multiIndex&>(fePtr->get_data_dofs()),MBTET,0,data.dataOnEntities[MBTET][0].getFieldData()); CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

// ** DoFS **

PetscErrorCode ForcesAndSurcesCore::getNodesFieldDofs(
  const string &field_name,FEDofMoFEMEntity_multiIndex &dofs,ublas::vector<const FEDofMoFEMEntity*> &nodes_field_dofs
) {
  PetscFunctionBegin;
  FEDofMoFEMEntity_multiIndex::index<Composite_Name_And_Type_mi_tag>::type::iterator dit,hi_dit;
  dit = dofs.get<Composite_Name_And_Type_mi_tag>().lower_bound(boost::make_tuple(field_name,MBVERTEX));
  hi_dit = dofs.get<Composite_Name_And_Type_mi_tag>().upper_bound(boost::make_tuple(field_name,MBVERTEX));
  nodes_field_dofs.resize(distance(dit,hi_dit));
  for(;dit!=hi_dit;dit++) {
    int side_number = dit->side_number_ptr->side_number;
    nodes_field_dofs[side_number*dit->get_max_rank()+dit->get_dof_rank()] = &*dit;
    int brother_side_number = dit->side_number_ptr->brother_side_number;
    if(brother_side_number!=-1) {
      if(nodes_field_dofs.size()<(unsigned int)(brother_side_number*dit->get_max_rank()+dit->get_max_rank())) {
        nodes_field_dofs.resize(brother_side_number*dit->get_max_rank()+dit->get_max_rank());
      }
      nodes_field_dofs[brother_side_number*dit->get_max_rank()+dit->get_dof_rank()] = &*dit;
    }
  }
  PetscFunctionReturn(0);
}

PetscErrorCode ForcesAndSurcesCore::getNodesFieldDofs(DataForcesAndSurcesCore &data,const string &field_name) {
    PetscFunctionBegin;
    ierr = getNodesFieldDofs(
      field_name,const_cast<FEDofMoFEMEntity_multiIndex&>(fePtr->get_data_dofs()),data.dataOnEntities[MBVERTEX][0].getFieldDofs()); CHKERRQ(ierr);
    PetscFunctionReturn(0);
  }

PetscErrorCode ForcesAndSurcesCore::getTypeFieldDofs(
  const string &field_name,FEDofMoFEMEntity_multiIndex &dofs,EntityType type,int side_number,ublas::vector<const FEDofMoFEMEntity*> &ent_field_dofs) {
    PetscFunctionBegin;
    FEDofMoFEMEntity_multiIndex::index<Composite_Name_Type_And_Side_Number_mi_tag>::type::iterator dit,hi_dit;
    dit = dofs.get<Composite_Name_Type_And_Side_Number_mi_tag>().lower_bound(boost::make_tuple(field_name,type,side_number));
    hi_dit = dofs.get<Composite_Name_Type_And_Side_Number_mi_tag>().upper_bound(boost::make_tuple(field_name,type,side_number));
    ent_field_dofs.resize(0);
    for(;dit!=hi_dit;dit++) {
      ent_field_dofs.resize(dit->get_nb_dofs_on_ent());
      ent_field_dofs[dit->get_EntDofIdx()] = &*dit;
    }
    PetscFunctionReturn(0);
  }

PetscErrorCode ForcesAndSurcesCore::getTypeFieldDofs(
  const string &field_name,FEDofMoFEMEntity_multiIndex &dofs,EntityType type,
  boost::ptr_vector<DataForcesAndSurcesCore::EntData> &data) {
  PetscFunctionBegin;
  SideNumber_multiIndex& side_table = const_cast<SideNumber_multiIndex&>(fePtr->get_side_number_table());
  if(data.size() < side_table.get<2>().count(type)) {
    SETERRQ(PETSC_COMM_SELF,MOFEM_DATA_INCONSISTENCT,"data inconsistency");
  }
  SideNumber_multiIndex::nth_index<2>::type::iterator siit = side_table.get<2>().lower_bound(type);
  SideNumber_multiIndex::nth_index<2>::type::iterator hi_siit = side_table.get<2>().upper_bound(type);
  for(;siit!=hi_siit;siit++) {
    ierr = getTypeFieldDofs(field_name,dofs,type,siit->side_number,data[siit->side_number].getFieldDofs()); CHKERRQ(ierr);
    if(siit->brother_side_number!=-1) {
      ierr = getTypeFieldDofs(field_name,dofs,type,siit->side_number,data[siit->brother_side_number].getFieldDofs()); CHKERRQ(ierr);
    }
  }
  PetscFunctionReturn(0);
}

PetscErrorCode ForcesAndSurcesCore::getNoFieldFieldDofs(
  const string &field_name,FEDofMoFEMEntity_multiIndex &dofs,ublas::vector<const FEDofMoFEMEntity*> &nodes_dofs
) {
  PetscFunctionBegin;
  FEDofMoFEMEntity_multiIndex::index<FieldName_mi_tag>::type::iterator dit,hi_dit;
  dit = dofs.get<FieldName_mi_tag>().lower_bound(field_name);
  hi_dit = dofs.get<FieldName_mi_tag>().upper_bound(field_name);
  nodes_dofs.resize(distance(dit,hi_dit));
  for(;dit!=hi_dit;dit++) {
    nodes_dofs[dit->get_dof_rank()] = &*dit;
  }
  PetscFunctionReturn(0);
}

PetscErrorCode ForcesAndSurcesCore::getNoFieldFieldDofs(
  DataForcesAndSurcesCore &data,const string &field_name
) {
  PetscFunctionBegin;
  if(data.dataOnEntities[MBENTITYSET].size() == 0) {
    SETERRQ(PETSC_COMM_SELF,MOFEM_DATA_INCONSISTENCT,"data inconsistency");
  }
  ierr = getNoFieldFieldDofs(
    field_name,const_cast<FEDofMoFEMEntity_multiIndex&>(fePtr->get_data_dofs()),data.dataOnEntities[MBENTITYSET][0].getFieldDofs()
  ); CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

PetscErrorCode ForcesAndSurcesCore::getEdgesFieldDofs(DataForcesAndSurcesCore &data,const string &field_name) {
  PetscFunctionBegin;
  ierr = getTypeFieldDofs(
    field_name,const_cast<FEDofMoFEMEntity_multiIndex&>(fePtr->get_data_dofs()),MBEDGE,data.dataOnEntities[MBEDGE]
  ); CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

PetscErrorCode ForcesAndSurcesCore::getTrisFieldDofs(DataForcesAndSurcesCore &data,const string &field_name) {
  PetscFunctionBegin;
  ierr = getTypeFieldDofs(
    field_name,const_cast<FEDofMoFEMEntity_multiIndex&>(fePtr->get_data_dofs()),MBTRI,data.dataOnEntities[MBTRI]
  ); CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

PetscErrorCode ForcesAndSurcesCore::getTetsFieldDofs(DataForcesAndSurcesCore &data,const string &field_name) {
  PetscFunctionBegin;
  if(data.dataOnEntities[MBTET].size() == 0) {
    SETERRQ(PETSC_COMM_SELF,MOFEM_DATA_INCONSISTENCT,"data inconsistency");
  }
  ierr = getTypeFieldDofs(
    field_name,const_cast<FEDofMoFEMEntity_multiIndex&>(fePtr->get_data_dofs()),MBTET,0,data.dataOnEntities[MBTET][0].getFieldDofs()
  ); CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

// ** Face **

PetscErrorCode ForcesAndSurcesCore::getFaceNodes(DataForcesAndSurcesCore &data) {
  PetscFunctionBegin;
  //PetscAttachDebugger();
  data.facesNodes.resize(4,3);
  SideNumber_multiIndex& side_table = const_cast<SideNumber_multiIndex&>(fePtr->get_side_number_table());
  const int cannonical_face_sense_p1[4][3] = { {0,1,3}, {1,2,3}, {0,3,2}/**/, {0,2,1}/**/ }; //secon index is offset (positive sense)
  const int cannonical_face_sense_m1[4][3] = { {0,3,1}, {1,3,2}, {0,2,3}, {0,1,2} }; //second index is offset (negative sense)
  SideNumber_multiIndex::nth_index<1>::type::iterator siit = side_table.get<1>().lower_bound(boost::make_tuple(MBTRI,0));
  SideNumber_multiIndex::nth_index<1>::type::iterator hi_siit = side_table.get<1>().upper_bound(boost::make_tuple(MBTRI,4));
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
      rval = mField.get_moab().get_connectivity(ent,conn_tet,num_nodes_tet,true); CHKERR_PETSC(rval);
      if(num_nodes_tet != 4) SETERRQ(PETSC_COMM_SELF,MOFEM_DATA_INCONSISTENCT,"data inconsistency");
      int num_nodes_face;
      const EntityHandle *conn_face;
      rval = mField.get_moab().get_connectivity(side->ent,conn_face,num_nodes_face,true); CHKERR_PETSC(rval);
      if(num_nodes_face != 3) SETERRQ(PETSC_COMM_SELF,1,"data inconsistency");
      if(conn_face[0] != conn_tet[data.facesNodes(side->side_number,0)]) SETERRQ(PETSC_COMM_SELF,MOFEM_DATA_INCONSISTENCT,"data inconsistency");
      if(conn_face[1] != conn_tet[data.facesNodes(side->side_number,1)]) SETERRQ(PETSC_COMM_SELF,MOFEM_DATA_INCONSISTENCT,"data inconsistency");
      if(conn_face[2] != conn_tet[data.facesNodes(side->side_number,2)]) SETERRQ(PETSC_COMM_SELF,MOFEM_DATA_INCONSISTENCT,"data inconsistency");
    }
  }
  PetscFunctionReturn(0);
}

// ** Space **

PetscErrorCode ForcesAndSurcesCore::getSpacesOnEntities(DataForcesAndSurcesCore &data) {
  PetscFunctionBegin;

  try {

  for(_IT_GET_FEDATA_DOFS_FOR_LOOP_(this,dof)) {
    data.spacesOnEntities[dof->get_ent_type()].set(dof->get_space());
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
    const double *G_X,const double *G_Y,const double *G_Z,const int G_DIM) {
    PetscFunctionBegin;

    data.dataOnEntities[MBVERTEX][0].getN().resize(G_DIM,4);
    ierr = ShapeMBTET(&*data.dataOnEntities[MBVERTEX][0].getN().data().begin(),G_X,G_Y,G_Z,G_DIM); CHKERRQ(ierr);
    data.dataOnEntities[MBVERTEX][0].getDiffN().resize(4,3);
    ierr = ShapeDiffMBTET(&*data.dataOnEntities[MBVERTEX][0].getDiffN().data().begin()); CHKERRQ(ierr);

    if((data.spacesOnEntities[MBEDGE]).test(H1)) {

    //edges
    if(data.dataOnEntities[MBEDGE].size()!=6) {
      SETERRQ(PETSC_COMM_SELF,MOFEM_DATA_INCONSISTENCT,"data inconsistency");
    }
    int _sense_[6],_order_[6];
    double *_H1edgeN_[6],*_diffH1edgeN_[6];
    for(int ee = 0;ee<6;ee++) {
      if(data.dataOnEntities[MBEDGE][ee].getSense() == 0) {
        SETERRQ(PETSC_COMM_SELF,MOFEM_DATA_INCONSISTENCT,"data inconsistency");
      }
      _sense_[ee] = data.dataOnEntities[MBEDGE][ee].getSense();
      _order_[ee] = data.dataOnEntities[MBEDGE][ee].getOrder();
      int nb_dofs = NBEDGE_H1(data.dataOnEntities[MBEDGE][ee].getOrder());
      data.dataOnEntities[MBEDGE][ee].getN().resize(G_DIM,nb_dofs);
      data.dataOnEntities[MBEDGE][ee].getDiffN().resize(G_DIM,3*nb_dofs);
      _H1edgeN_[ee] = &*data.dataOnEntities[MBEDGE][ee].getN().data().begin();
      _diffH1edgeN_[ee] = &*data.dataOnEntities[MBEDGE][ee].getDiffN().data().begin();
    }
    ierr = H1_EdgeShapeFunctions_MBTET(
      _sense_,_order_,
      &*data.dataOnEntities[MBVERTEX][0].getN().data().begin(),&*data.dataOnEntities[MBVERTEX][0].getDiffN().data().begin(),
      _H1edgeN_,_diffH1edgeN_,G_DIM); CHKERRQ(ierr);

    //faces
    if(data.dataOnEntities[MBTRI].size()!=4) {
      SETERRQ(PETSC_COMM_SELF,MOFEM_DATA_INCONSISTENCT,"data inconsistency");
    }
    double *_H1faceN_[4],*_diffH1faceN_[4];
    for(int ff = 0;ff<4;ff++) {
      if(data.dataOnEntities[MBTRI][ff].getSense() == 0) {
        SETERRQ(PETSC_COMM_SELF,MOFEM_DATA_INCONSISTENCT,"data inconsistency");
      }
      int nb_dofs = NBFACE_H1(data.dataOnEntities[MBTRI][ff].getOrder());
      _order_[ff] = data.dataOnEntities[MBTRI][ff].getOrder();
      data.dataOnEntities[MBTRI][ff].getN().resize(G_DIM,nb_dofs);
      data.dataOnEntities[MBTRI][ff].getDiffN().resize(G_DIM,3*nb_dofs);
      _H1faceN_[ff] = &*data.dataOnEntities[MBTRI][ff].getN().data().begin();
      _diffH1faceN_[ff] = &*data.dataOnEntities[MBTRI][ff].getDiffN().data().begin();
    }
    if(data.facesNodes.size1() != 4) {
      SETERRQ(PETSC_COMM_SELF,MOFEM_DATA_INCONSISTENCT,"data inconsistency");
    }
    if(data.facesNodes.size2() != 3) {
      SETERRQ(PETSC_COMM_SELF,MOFEM_DATA_INCONSISTENCT,"data inconsistency");
    }
    ierr = H1_FaceShapeFunctions_MBTET(
      &*data.facesNodes.data().begin(),_order_,
      &*data.dataOnEntities[MBVERTEX][0].getN().data().begin(),&*data.dataOnEntities[MBVERTEX][0].getDiffN().data().begin(),
      _H1faceN_,_diffH1faceN_,G_DIM); CHKERRQ(ierr);

    //volume
    int nb_vol_dofs = NBVOLUME_H1(data.dataOnEntities[MBTET][0].getOrder());
    data.dataOnEntities[MBTET][0].getN().resize(G_DIM,nb_vol_dofs);
    data.dataOnEntities[MBTET][0].getDiffN().resize(G_DIM,3*nb_vol_dofs);
    ierr = H1_VolumeShapeFunctions_MBTET(
      data.dataOnEntities[MBTET][0].getOrder(),&*data.dataOnEntities[MBVERTEX][0].getN().data().begin(),&*data.dataOnEntities[MBVERTEX][0].getDiffN().data().begin(),
      &*data.dataOnEntities[MBTET][0].getN().data().begin(),&*data.dataOnEntities[MBTET][0].getDiffN().data().begin(),G_DIM
    ); CHKERRQ(ierr);

    }

    PetscFunctionReturn(0);
  }

PetscErrorCode ForcesAndSurcesCore::shapeTETFunctions_L2(
    DataForcesAndSurcesCore &data,
    const double *G_X,const double *G_Y,const double *G_Z,const int G_DIM) {
    PetscFunctionBegin;

    data.dataOnEntities[MBVERTEX][0].getN().resize(G_DIM,4);
    ierr = ShapeMBTET(&*data.dataOnEntities[MBVERTEX][0].getN().data().begin(),G_X,G_Y,G_Z,G_DIM); CHKERRQ(ierr);
    data.dataOnEntities[MBVERTEX][0].getDiffN().resize(4,3);
    ierr = ShapeDiffMBTET(&*data.dataOnEntities[MBVERTEX][0].getDiffN().data().begin()); CHKERRQ(ierr);

    data.dataOnEntities[MBTET][0].getN().resize(G_DIM,NBVOLUME_L2(data.dataOnEntities[MBTET][0].getOrder()));
    data.dataOnEntities[MBTET][0].getDiffN().resize(G_DIM,3*NBVOLUME_L2(data.dataOnEntities[MBTET][0].getOrder()));

    ierr = L2_ShapeFunctions_MBTET(
      data.dataOnEntities[MBTET][0].getOrder(),
      &*data.dataOnEntities[MBVERTEX][0].getN().data().begin(),&*data.dataOnEntities[MBVERTEX][0].getDiffN().data().begin(),
      &*data.dataOnEntities[MBTET][0].getN().data().begin(),
      &*data.dataOnEntities[MBTET][0].getDiffN().data().begin(),
      G_DIM); CHKERRQ(ierr);

    PetscFunctionReturn(0);
  }

PetscErrorCode ForcesAndSurcesCore::shapeTETFunctions_Hdiv(
    DataForcesAndSurcesCore &data,
    const double *G_X,const double *G_Y,const double *G_Z,const int G_DIM) {
    PetscFunctionBegin;

    //calculate shape function for tet, needed to construct shape functions for h_div

    if(data.dataOnEntities[MBVERTEX][0].getN().size1() != (unsigned int)G_DIM) {
      data.dataOnEntities[MBVERTEX][0].getN().resize(G_DIM,4);
      ierr = ShapeMBTET(&*data.dataOnEntities[MBVERTEX][0].getN().data().begin(),G_X,G_Y,G_Z,G_DIM); CHKERRQ(ierr);
    }
    //that is cheep to calate, no harm done if recalculated
    data.dataOnEntities[MBVERTEX][0].getDiffN().resize(4,3);
    ierr = ShapeDiffMBTET(&*data.dataOnEntities[MBVERTEX][0].getDiffN().data().begin()); CHKERRQ(ierr);

    //face shape functions

    double *phi_f_e[4][3];
    double *phi_f[4];
    double *diff_phi_f_e[4][3];
    double *diff_phi_f[4];

    N_face_edge.resize(4,3);
    N_face_bubble.resize(4);
    diffN_face_edge.resize(4,3);
    diffN_face_bubble.resize(4);

    int faces_order[4];
    for(int ff = 0;ff<4;ff++) {
      if(data.dataOnEntities[MBTRI][ff].getSense() == 0) {
	SETERRQ(PETSC_COMM_SELF,MOFEM_DATA_INCONSISTENCT,"data inconsistency");
      }
      faces_order[ff] = data.dataOnEntities[MBTRI][ff].getOrder();
      //three edges on face
      for(int ee = 0;ee<3;ee++) {
        N_face_edge(ff,ee).resize(G_DIM,3*NBFACE_EDGE_HDIV(faces_order[ff]));
        diffN_face_edge(ff,ee).resize(G_DIM,9*NBFACE_EDGE_HDIV(faces_order[ff]));
        phi_f_e[ff][ee] = &((N_face_edge(ff,ee))(0,0));
        diff_phi_f_e[ff][ee] = &((diffN_face_edge(ff,ee))(0,0));
      }
      N_face_bubble[ff].resize(G_DIM,3*NBFACE_FACE_HDIV(faces_order[ff]));
      diffN_face_bubble[ff].resize(G_DIM,9*NBFACE_FACE_HDIV(faces_order[ff]));
      phi_f[ff] = &*(N_face_bubble[ff].data().begin());
      diff_phi_f[ff] = &*(diffN_face_bubble[ff].data().begin());
    }

    ierr = Hdiv_EdgeFaceShapeFunctions_MBTET(
      &data.facesNodes(0,0),faces_order,
      &data.dataOnEntities[MBVERTEX][0].getN()(0,0),
      &data.dataOnEntities[MBVERTEX][0].getDiffN()(0,0),
      phi_f_e,diff_phi_f_e,G_DIM); CHKERRQ(ierr);

    ierr = Hdiv_FaceBubbleShapeFunctions_MBTET(
      &data.facesNodes(0,0),faces_order,
      &data.dataOnEntities[MBVERTEX][0].getN()(0,0),
      &data.dataOnEntities[MBVERTEX][0].getDiffN()(0,0),
      phi_f,diff_phi_f,G_DIM); CHKERRQ(ierr);

    //volume shape functions

    double *phi_v_e[6];
    double *phi_v_f[4];
    double *phi_v;
    double *diff_phi_v_e[6];
    double *diff_phi_v_f[4];
    double *diff_phi_v;

    int volume_order = data.dataOnEntities[MBTET][0].getOrder();
    double coords[] = { 0,0,0, 1,0,0, 0,1,0, 0,0,1 };

    N_volume_edge.resize(6);
    diffN_volume_edge.resize(6);
    for(int ee = 0;ee<6;ee++) {
      N_volume_edge[ee].resize(G_DIM,3*NBVOLUME_EDGE_HDIV(volume_order));
      diffN_volume_edge[ee].resize(G_DIM,9*NBVOLUME_EDGE_HDIV(volume_order));
      phi_v_e[ee] = &*(N_volume_edge[ee].data().begin());
      diff_phi_v_e[ee] = &*(diffN_volume_edge[ee].data().begin());
    }
    ierr = Hdiv_EdgeBasedVolumeShapeFunctions_MBTET(
      volume_order,coords,&data.dataOnEntities[MBVERTEX][0].getN()(0,0),
      &data.dataOnEntities[MBVERTEX][0].getDiffN()(0,0),
      phi_v_e,diff_phi_v_e,G_DIM); CHKERRQ(ierr);

    N_volume_face.resize(4);
    diffN_volume_face.resize(4);
    for(int ff = 0;ff<4;ff++) {
      N_volume_face[ff].resize(G_DIM,3*NBVOLUME_FACE_HDIV(volume_order));
      diffN_volume_face[ff].resize(G_DIM,9*NBVOLUME_FACE_HDIV(volume_order));
      phi_v_f[ff] = &*(N_volume_face[ff].data().begin());
      diff_phi_v_f[ff] = &*(diffN_volume_face[ff].data().begin());
    }
    ierr = Hdiv_FaceBasedVolumeShapeFunctions_MBTET(
      volume_order,coords,&data.dataOnEntities[MBVERTEX][0].getN()(0,0),
      &data.dataOnEntities[MBVERTEX][0].getDiffN()(0,0),
      phi_v_f,diff_phi_v_f,G_DIM); CHKERRQ(ierr);

    N_volume_bubble.resize(G_DIM,3*NBVOLUME_VOLUME_HDIV(volume_order));
    diffN_volume_bubble.resize(G_DIM,9*NBVOLUME_VOLUME_HDIV(volume_order));
    phi_v = &*(N_volume_bubble.data().begin());
    diff_phi_v = &*(diffN_volume_bubble.data().begin());
    ierr = Hdiv_VolumeBubbleShapeFunctions_MBTET(
      volume_order,coords,&data.dataOnEntities[MBVERTEX][0].getN()(0,0),
      &data.dataOnEntities[MBVERTEX][0].getDiffN()(0,0),
      phi_v,diff_phi_v,G_DIM); CHKERRQ(ierr);

    // Set shape functions into data strucrure Shape functions hast to be put
    // in arrays in order which guarantee hierarhical series of digrees of
    // freedom, i.e. in other words dofs form sub-entities has to be group
    // by order.

    //faces
    if(data.dataOnEntities[MBTRI].size()!=4) {
      SETERRQ(PETSC_COMM_SELF,MOFEM_DATA_INCONSISTENCT,"data inconsistency");
    }
    for(int ff = 0;ff<4;ff++) {
      data.dataOnEntities[MBTRI][ff].getHdivN().resize(G_DIM,3*NBFACE_HDIV(faces_order[ff]),0);
      data.dataOnEntities[MBTRI][ff].getDiffHdivN().resize(G_DIM,9*NBFACE_HDIV(faces_order[ff]),0);
      int col = 0,diff_col = 0;
      for(int oo = 0;oo<faces_order[ff];oo++) {
        for(int ee = 0;ee<3;ee++) {
          //values
          for(int dd = 3*NBFACE_EDGE_HDIV(oo);dd<3*NBFACE_EDGE_HDIV(oo+1);dd++,col++) {
            for(int gg = 0;gg<G_DIM;gg++) {
              data.dataOnEntities[MBTRI][ff].getHdivN()(gg,col) = N_face_edge(ff,ee)(gg,dd);
            }
          }
          //direvatives
          for(int dd = 9*NBFACE_EDGE_HDIV(oo);dd<9*NBFACE_EDGE_HDIV(oo+1);dd++,diff_col++) {
            for(int gg = 0;gg<G_DIM;gg++) {
              data.dataOnEntities[MBTRI][ff].getDiffHdivN()(gg,diff_col) = diffN_face_edge(ff,ee)(gg,dd);
            }
          }
        }
        //values
        for(int dd = 3*NBFACE_FACE_HDIV(oo);dd<3*NBFACE_FACE_HDIV(oo+1);dd++,col++) {
          for(int gg = 0;gg<G_DIM;gg++) {
            data.dataOnEntities[MBTRI][ff].getHdivN()(gg,col) = N_face_bubble[ff](gg,dd);
          }
        }
        //direvatives
        for(int dd = 9*NBFACE_FACE_HDIV(oo);dd<9*NBFACE_FACE_HDIV(oo+1);dd++,diff_col++) {
          for(int gg = 0;gg<G_DIM;gg++) {
            data.dataOnEntities[MBTRI][ff].getDiffHdivN()(gg,diff_col) = diffN_face_bubble[ff](gg,dd);
          }
        }
      }
    }

    //volume
    int col = 0,diff_col = 0;
    data.dataOnEntities[MBTET][0].getHdivN().resize(G_DIM,3*NBVOLUME_HDIV(volume_order),0);
    data.dataOnEntities[MBTET][0].getDiffHdivN().resize(G_DIM,9*NBVOLUME_HDIV(volume_order),0);
    for(int oo = 0;oo<volume_order;oo++) {
      for(int ee = 0;ee<6;ee++) {
	//values
	for(int dd = 3*NBVOLUME_EDGE_HDIV(oo);dd<3*NBVOLUME_EDGE_HDIV(oo+1);dd++,col++) {
	  for(int gg = 0;gg<G_DIM;gg++) {
	    data.dataOnEntities[MBTET][0].getHdivN()(gg,col) = N_volume_edge[ee](gg,dd);
	  }
	}
	//direvatives
	for(int dd = 9*NBVOLUME_EDGE_HDIV(oo);dd<9*NBVOLUME_EDGE_HDIV(oo+1);dd++,diff_col++) {
	  for(int gg = 0;gg<G_DIM;gg++) {
	    data.dataOnEntities[MBTET][0].getDiffHdivN()(gg,diff_col) = diffN_volume_edge[ee](gg,dd);
	  }
	}
      }
      for(int ff = 0;ff<4;ff++) {
	//values
	for(int dd = 3*NBVOLUME_FACE_HDIV(oo);dd<3*NBVOLUME_FACE_HDIV(oo+1);dd++,col++) {
	  for(int gg = 0;gg<G_DIM;gg++) {
	    data.dataOnEntities[MBTET][0].getHdivN()(gg,col) = N_volume_face[ff](gg,dd);
	  }
	}
	//direvatives
	for(int dd = 9*NBVOLUME_FACE_HDIV(oo);dd<9*NBVOLUME_FACE_HDIV(oo+1);dd++,diff_col++) {
	  for(int gg = 0;gg<G_DIM;gg++) {
	    data.dataOnEntities[MBTET][0].getDiffHdivN()(gg,diff_col) = diffN_volume_face[ff](gg,dd);
	  }
	}
      }
      //values
      for(int dd = 3*NBVOLUME_VOLUME_HDIV(oo);dd<3*NBVOLUME_VOLUME_HDIV(oo+1);dd++,col++) {
	for(int gg = 0;gg<G_DIM;gg++) {
	  data.dataOnEntities[MBTET][0].getHdivN()(gg,col) = N_volume_bubble(gg,dd);
	}
      }
      //direvatives
      for(int dd = 9*NBVOLUME_VOLUME_HDIV(oo);dd<9*NBVOLUME_VOLUME_HDIV(oo+1);dd++,diff_col++) {
	for(int gg = 0;gg<G_DIM;gg++) {
	  data.dataOnEntities[MBTET][0].getDiffHdivN()(gg,diff_col) = diffN_volume_bubble(gg,dd);
	}
      }

    }

    PetscFunctionReturn(0);
  }

PetscErrorCode ForcesAndSurcesCore::shapeTRIFunctions_H1(
  DataForcesAndSurcesCore &data,
  const double *G_X,const double *G_Y,const int G_DIM) {
  PetscFunctionBegin;

  data.dataOnEntities[MBVERTEX][0].getN().resize(G_DIM,3);
  ierr = ShapeMBTRI(&*data.dataOnEntities[MBVERTEX][0].getN().data().begin(),G_X,G_Y,G_DIM); CHKERRQ(ierr);
  data.dataOnEntities[MBVERTEX][0].getDiffN().resize(3,2);
  ierr = ShapeDiffMBTRI(&*data.dataOnEntities[MBVERTEX][0].getDiffN().data().begin()); CHKERRQ(ierr);

  if((data.spacesOnEntities[MBEDGE]).test(H1)) {

  //edges
  if(data.dataOnEntities[MBEDGE].size()!=3) {
    SETERRQ(PETSC_COMM_SELF,MOFEM_DATA_INCONSISTENCT,"data inconsistency");
  }
  int _sense_[3],_order_[3];
  double *_H1edgeN_[3],*_diffH1edgeN_[3];
  for(int ee = 0;ee<3;ee++) {
    if(data.dataOnEntities[MBEDGE][ee].getSense() == 0) {
	SETERRQ(PETSC_COMM_SELF,MOFEM_DATA_INCONSISTENCT,"data inconsistency");
    }
    _sense_[ee] = data.dataOnEntities[MBEDGE][ee].getSense();
    _order_[ee] = data.dataOnEntities[MBEDGE][ee].getOrder();
    int nb_dofs = NBEDGE_H1(data.dataOnEntities[MBEDGE][ee].getOrder());
    data.dataOnEntities[MBEDGE][ee].getN().resize(G_DIM,nb_dofs);
    data.dataOnEntities[MBEDGE][ee].getDiffN().resize(G_DIM,2*nb_dofs);
    _H1edgeN_[ee] = &*data.dataOnEntities[MBEDGE][ee].getN().data().begin();
    _diffH1edgeN_[ee] = &*data.dataOnEntities[MBEDGE][ee].getDiffN().data().begin();
  }
  ierr = H1_EdgeShapeFunctions_MBTRI(_sense_,_order_,
    &*data.dataOnEntities[MBVERTEX][0].getN().data().begin(),&*data.dataOnEntities[MBVERTEX][0].getDiffN().data().begin(),
    _H1edgeN_,_diffH1edgeN_,G_DIM); CHKERRQ(ierr);

  //face
  if(data.dataOnEntities[MBTRI].size()!=1) {
    SETERRQ(PETSC_COMM_SELF,MOFEM_DATA_INCONSISTENCT,"data inconsistency");
  }
  int nb_dofs = NBFACE_H1(data.dataOnEntities[MBTRI][0].getOrder());
  data.dataOnEntities[MBTRI][0].getN().resize(G_DIM,nb_dofs);
  data.dataOnEntities[MBTRI][0].getDiffN().resize(G_DIM,2*nb_dofs);
  const int face_nodes[] = { 0,1,2 };
  ierr = H1_FaceShapeFunctions_MBTRI(
    face_nodes,data.dataOnEntities[MBTRI][0].getOrder(),
    &*data.dataOnEntities[MBVERTEX][0].getN().data().begin(),&*data.dataOnEntities[MBVERTEX][0].getDiffN().data().begin(),
    &*data.dataOnEntities[MBTRI][0].getN().data().begin(),&*data.dataOnEntities[MBTRI][0].getDiffN().data().begin(),
    G_DIM); CHKERRQ(ierr);

  }

  PetscFunctionReturn(0);
}

PetscErrorCode ForcesAndSurcesCore::shapeTRIFunctions_Hdiv(
    DataForcesAndSurcesCore &data,
    const double *G_X,const double *G_Y,const int G_DIM) {
  PetscFunctionBegin;


  data.dataOnEntities[MBVERTEX][0].getN().resize(G_DIM,3);
  ierr = ShapeMBTRI(&*data.dataOnEntities[MBVERTEX][0].getN().data().begin(),G_X,G_Y,G_DIM); CHKERRQ(ierr);

  double *PHI_f_e[3];
  double *PHI_f;

  N_face_edge.resize(1,3);
  N_face_bubble.resize(1);
  int face_order = data.dataOnEntities[MBTRI][0].getOrder();
  //three edges on face
  for(int ee = 0;ee<3;ee++) {
    N_face_edge(0,ee).resize(G_DIM,3*NBFACE_EDGE_HDIV(face_order));
    PHI_f_e[ee] = &((N_face_edge(0,ee))(0,0));
  }
  N_face_bubble[0].resize(G_DIM,3*NBFACE_FACE_HDIV(face_order));
  PHI_f = &*(N_face_bubble[0].data().begin());

  int face_nodes[3] = { 0,1,2 };
  ierr = Hdiv_EdgeFaceShapeFunctions_MBTET_ON_FACE(face_nodes,face_order,
    &data.dataOnEntities[MBVERTEX][0].getN()(0,0),NULL,
    PHI_f_e,NULL,G_DIM,3); CHKERRQ(ierr);
  ierr = Hdiv_FaceBubbleShapeFunctions_MBTET_ON_FACE(face_nodes,face_order,
    &data.dataOnEntities[MBVERTEX][0].getN()(0,0),NULL,
    PHI_f,NULL,G_DIM,3); CHKERRQ(ierr);

  // set shape functions into data strucrure

  if(data.dataOnEntities[MBTRI].size()!=1) {
    SETERRQ(PETSC_COMM_SELF,MOFEM_DATA_INCONSISTENCT,"data inconsistency");
  }
  data.dataOnEntities[MBTRI][0].getHdivN().resize(G_DIM,3*NBFACE_HDIV(face_order),0);
  int col = 0;
  for(int oo = 0;oo<face_order;oo++) {
    for(int ee = 0;ee<3;ee++) {
      for(int dd = 3*NBFACE_EDGE_HDIV(oo);dd<3*NBFACE_EDGE_HDIV(oo+1);dd++,col++) {
	for(int gg = 0;gg<G_DIM;gg++) {
	  data.dataOnEntities[MBTRI][0].getHdivN()(gg,col) = N_face_edge(0,ee)(gg,dd);
	}
      }
    }
    for(int dd = 3*NBFACE_FACE_HDIV(oo);dd<3*NBFACE_FACE_HDIV(oo+1);dd++,col++) {
      for(int gg = 0;gg<G_DIM;gg++) {
	data.dataOnEntities[MBTRI][0].getHdivN()(gg,col) = N_face_bubble[0](gg,dd);
      }
    }
  }

  PetscFunctionReturn(0);
}

PetscErrorCode ForcesAndSurcesCore::shapeEDGEFunctions_H1(DataForcesAndSurcesCore &data,const double *G_X,const int G_DIM) {
  PetscFunctionBegin;

  data.dataOnEntities[MBVERTEX][0].getN().resize(G_DIM,2);
  ierr = ShapeMBEDGE(&*data.dataOnEntities[MBVERTEX][0].getN().data().begin(),G_X,G_DIM); CHKERRQ(ierr);
  data.dataOnEntities[MBVERTEX][0].getDiffN().resize(2,1);
  ierr = ShapeDiffMBEDGE(&*data.dataOnEntities[MBVERTEX][0].getDiffN().data().begin()); CHKERRQ(ierr);

  //cerr << data.dataOnEntities[MBVERTEX][0].getN() << endl;
  //cerr << data.dataOnEntities[MBVERTEX][0].getDiffN() << endl;

  int order = data.dataOnEntities[MBEDGE][0].getOrder();
  data.dataOnEntities[MBEDGE][0].getN().resize(G_DIM,NBEDGE_H1(order));
  data.dataOnEntities[MBEDGE][0].getDiffN().resize(G_DIM,NBEDGE_H1(order));
  if(data.dataOnEntities[MBEDGE][0].getOrder()>1) {
    double diff_s = 0.5;
    for(int gg = 0;gg<G_DIM;gg++) {
      double s = 2*G_X[gg]-1;
      ierr = Lagrange_basis(NBEDGE_H1(order)-1,s,&diff_s,
	&data.dataOnEntities[MBEDGE][0].getN()(gg,0),&data.dataOnEntities[MBEDGE][0].getDiffN()(gg,0),1); CHKERRQ(ierr);
      for(unsigned int pp = 0;pp<data.dataOnEntities[MBEDGE][0].getN().size2();pp++) {
	double L = data.dataOnEntities[MBEDGE][0].getN()(gg,pp);
	double diffL = data.dataOnEntities[MBEDGE][0].getDiffN()(gg,pp);
	data.dataOnEntities[MBEDGE][0].getN()(gg,pp) = data.dataOnEntities[MBVERTEX][0].getN()(gg,0)*data.dataOnEntities[MBVERTEX][0].getN()(gg,1)*L;
	data.dataOnEntities[MBEDGE][0].getDiffN()(gg,pp) =
	  ((+1.)*data.dataOnEntities[MBVERTEX][0].getN()(gg,1)+data.dataOnEntities[MBVERTEX][0].getN()(gg,0)*(-1.))*L + data.dataOnEntities[MBVERTEX][0].getN()(gg,0)*data.dataOnEntities[MBVERTEX][0].getN()(gg,1)*diffL;
      }
    }
  }

  //cerr << data.dataOnEntities[MBEDGE][0].getN() << endl;
  //cerr << data.dataOnEntities[MBEDGE][0].getDiffN() << endl;

  PetscFunctionReturn(0);
}

PetscErrorCode ForcesAndSurcesCore::shapeFlatPRISMFunctions_H1(
  DataForcesAndSurcesCore &data,
  const double *G_X,const double *G_Y,const int G_DIM) {
  PetscFunctionBegin;

  try {

  SideNumber_multiIndex& side_table = const_cast<SideNumber_multiIndex&>(fePtr->get_side_number_table());
  SideNumber_multiIndex::nth_index<1>::type::iterator siit3 = side_table.get<1>().find(boost::make_tuple(MBTRI,3));
  SideNumber_multiIndex::nth_index<1>::type::iterator siit4 = side_table.get<1>().find(boost::make_tuple(MBTRI,4));
  if(siit3==side_table.get<1>().end()) SETERRQ(PETSC_COMM_SELF,MOFEM_DATA_INCONSISTENCT,"data inconsistency");
  if(siit4==side_table.get<1>().end()) SETERRQ(PETSC_COMM_SELF,MOFEM_DATA_INCONSISTENCT,"data inconsistency");
  int num_nodes;
  const EntityHandle *conn_face3;
  const EntityHandle *conn_face4;
  rval = mField.get_moab().get_connectivity(siit3->ent,conn_face3,num_nodes,true); CHKERR_PETSC(rval);
  rval = mField.get_moab().get_connectivity(siit4->ent,conn_face4,num_nodes,true); CHKERR_PETSC(rval);

  ublas::matrix<double> N(G_DIM,3),diffN(3,2);
  ierr = ShapeMBTRI(&*N.data().begin(),G_X,G_Y,G_DIM); CHKERRQ(ierr);
  ierr = ShapeDiffMBTRI(&*diffN.data().begin()); CHKERRQ(ierr);
  data.dataOnEntities[MBVERTEX][0].getN().resize(G_DIM,6);
  data.dataOnEntities[MBVERTEX][0].getDiffN().resize(G_DIM,6*2);

  int face_nodes[2][3];
  for(int nn = 0;nn<3;nn++) {
    int side_number3 = fePtr->get_side_number_ptr(mField.get_moab(),conn_face3[nn])->side_number;
    int side_number4 = fePtr->get_side_number_ptr(mField.get_moab(),conn_face4[nn])->side_number;
    if(side_number3>2) {
      SETERRQ(PETSC_COMM_SELF,MOFEM_DATA_INCONSISTENCT,"data inconsistency");
    }
    if(side_number4 < 3) {
      side_number4 = fePtr->get_side_number_ptr(mField.get_moab(),conn_face4[nn])->brother_side_number;
      if(side_number4 == -1) {
	SETERRQ(PETSC_COMM_SELF,MOFEM_DATA_INCONSISTENCT,"data inconsistency");
      }
    }
    if(side_number4<3) {
      SETERRQ(PETSC_COMM_SELF,MOFEM_DATA_INCONSISTENCT,"data inconsistency");
    }
    face_nodes[0][nn] = side_number3;
    face_nodes[1][nn] = side_number4-3;
    for(int gg = 0;gg<G_DIM;gg++) {
      double val3 = N(gg,side_number3);
      double val3_x = diffN(side_number3,0);
      double val3_y = diffN(side_number3,1);
      data.dataOnEntities[MBVERTEX][0].getN()(gg,side_number3) = val3;
      data.dataOnEntities[MBVERTEX][0].getDiffN()(gg,2*side_number3+0) = val3_x;
      data.dataOnEntities[MBVERTEX][0].getDiffN()(gg,2*side_number3+1) = val3_y;
      double val4 = N(gg,side_number4-3);
      double val4_x = diffN(side_number4-3,0);
      double val4_y = diffN(side_number4-3,1);
      data.dataOnEntities[MBVERTEX][0].getN()(gg,side_number4) = val4;
      data.dataOnEntities[MBVERTEX][0].getDiffN()(gg,2*side_number4+0) = val4_x;
      data.dataOnEntities[MBVERTEX][0].getDiffN()(gg,2*side_number4+1) = val4_y;
    }
  }

  //edges
  if(data.dataOnEntities[MBEDGE].size()!=9) {
    SETERRQ(PETSC_COMM_SELF,MOFEM_DATA_INCONSISTENCT,"data inconsistency");
  }

  int valid_edges[] = { 1,1,1, 0,0,0, 1,1,1 };
  int sense[9],order[9];
  double *H1edgeN[9],*diffH1edgeN[9];
  for(int ee = 0;ee<9;ee++) {
    if(!valid_edges[ee]) continue;
    if(data.dataOnEntities[MBEDGE][ee].getSense() == 0) {
	SETERRQ(PETSC_COMM_SELF,MOFEM_DATA_INCONSISTENCT,"data inconsistency");
    }
    sense[ee] = data.dataOnEntities[MBEDGE][ee].getSense();
    order[ee] = data.dataOnEntities[MBEDGE][ee].getOrder();
    int nb_dofs = NBEDGE_H1(data.dataOnEntities[MBEDGE][ee].getOrder());
    data.dataOnEntities[MBEDGE][ee].getN().resize(G_DIM,nb_dofs);
    data.dataOnEntities[MBEDGE][ee].getDiffN().resize(G_DIM,2*nb_dofs);
    H1edgeN[ee] = &*data.dataOnEntities[MBEDGE][ee].getN().data().begin();
    diffH1edgeN[ee] = &*data.dataOnEntities[MBEDGE][ee].getDiffN().data().begin();
  }
  //shape functions on face 3
  ierr = H1_EdgeShapeFunctions_MBTRI(&sense[0],&order[0],
    &*N.data().begin(),&*diffN.data().begin(),
    &H1edgeN[0],&diffH1edgeN[0],G_DIM); CHKERRQ(ierr);
  //shape functions on face 4
  ierr = H1_EdgeShapeFunctions_MBTRI(&sense[6],&order[6],
    &*N.data().begin(),&*diffN.data().begin(),
    &H1edgeN[6],&diffH1edgeN[6],G_DIM); CHKERRQ(ierr);

  //face
  if(data.dataOnEntities[MBTRI].size()!=5) {
    SETERRQ(PETSC_COMM_SELF,MOFEM_DATA_INCONSISTENCT,"data inconsistency");
  }
  for(int ff = 3;ff<5;ff++) {
    int nb_dofs = NBFACE_H1(data.dataOnEntities[MBTRI][ff].getOrder());
    data.dataOnEntities[MBTRI][ff].getN().resize(G_DIM,nb_dofs);
    data.dataOnEntities[MBTRI][ff].getDiffN().resize(G_DIM,2*nb_dofs);
    ierr = H1_FaceShapeFunctions_MBTRI(
      face_nodes[ff-3],data.dataOnEntities[MBTRI][ff].getOrder(),
      &*N.data().begin(),&*diffN.data().begin(),
      &*data.dataOnEntities[MBTRI][ff].getN().data().begin(),&*data.dataOnEntities[MBTRI][ff].getDiffN().data().begin(),
      G_DIM); CHKERRQ(ierr);
  }

  } catch (exception& ex) {
    ostringstream ss;
    ss << "thorw in method: " << ex.what() << " at line " << __LINE__ << " in file " << __FILE__;
    SETERRQ(PETSC_COMM_SELF,MOFEM_STD_EXCEPTION_THROW,ss.str().c_str());
  }

  PetscFunctionReturn(0);
}

PetscErrorCode ForcesAndSurcesCore::shapeFlatPRISMFunctions_Hdiv(
    DataForcesAndSurcesCore &data,
    const double *G_X,const double *G_Y,const int G_DIM) {
  PetscFunctionBegin;
  SETERRQ(PETSC_COMM_SELF,MOFEM_NOT_IMPLEMENTED,"not yet implemented, i.e. flat prism with Hdiv space");
  PetscFunctionReturn(0);
}

// **** Data Operator ****

static PetscErrorCode get_porblem_row_indices(
  const ForcesAndSurcesCore *fe_ptr,const EntityType type,const int side,const string field_name,ublas::vector<int>& indices) {
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
  const ForcesAndSurcesCore *fe_ptr,const EntityType type,const int side,const string field_name,ublas::vector<int>& indices) {
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
  const string field_name,const EntityType type,const int side,ublas::vector<int>& indices) {
  PetscFunctionBegin;

  PetscErrorCode ierr;

  if(ptrFE == NULL) {
    SETERRQ(PETSC_COMM_SELF,MOFEM_DATA_INCONSISTENCT,"data inconsistency");
  }
  ierr = get_porblem_row_indices(ptrFE,type,side,field_name,indices); CHKERRQ(ierr);

  PetscFunctionReturn(0);
}

PetscErrorCode ForcesAndSurcesCore::UserDataOperator::getPorblemColIndices(
  const string field_name,const EntityType type,const int side,ublas::vector<int>& indices) {
  PetscFunctionBegin;

  PetscErrorCode ierr;

  if(ptrFE == NULL) {
    SETERRQ(PETSC_COMM_SELF,MOFEM_DATA_INCONSISTENCT,"data inconsistency");
  }
  ierr = get_porblem_col_indices(ptrFE,type,side,field_name,indices); CHKERRQ(ierr);

  PetscFunctionReturn(0);
}

// **** Tetrahedral ****

PetscErrorCode VolumeElementForcesAndSourcesCore::operator()() {
  PetscFunctionBegin;

  try {

    if(fePtr->get_ent_type() != MBTET) PetscFunctionReturn(0);

    ierr = getSpacesOnEntities(dataH1); CHKERRQ(ierr);

    //H1
    if((dataH1.spacesOnEntities[MBEDGE]).test(H1)) {
      ierr = getEdgesSense(dataH1); CHKERRQ(ierr);
      ierr = getTrisSense(dataH1); CHKERRQ(ierr);
      ierr = getEdgesOrder(dataH1,H1); CHKERRQ(ierr);
      ierr = getTrisOrder(dataH1,H1); CHKERRQ(ierr);
      ierr = getTetsOrder(dataH1,H1); CHKERRQ(ierr);
      ierr = getFaceNodes(dataH1); CHKERRQ(ierr);
    }

    //Hdiv
    if((dataH1.spacesOnEntities[MBTRI]).test(HDIV)) {
      ierr = getTrisSense(dataHdiv); CHKERRQ(ierr);
      ierr = getTrisOrder(dataHdiv,HDIV); CHKERRQ(ierr);
      ierr = getTetsOrder(dataHdiv,HDIV); CHKERRQ(ierr);
      ierr = getFaceNodes(dataHdiv); CHKERRQ(ierr);
    }

    //L2
    if((dataH1.spacesOnEntities[MBTET]).test(L2)) {
      ierr = getTetsOrder(dataL2,L2); CHKERRQ(ierr);
    }

    int order = 1;
    for(unsigned int ee = 0;ee<dataH1.dataOnEntities[MBEDGE].size();ee++) {
      order = max(order,dataH1.dataOnEntities[MBEDGE][ee].getOrder());
    }
    for(unsigned int ff = 0;ff<dataHdiv.dataOnEntities[MBTRI].size();ff++) {
      order = max(order,dataHdiv.dataOnEntities[MBTRI][ff].getOrder());
    }
    order = max(order,dataL2.dataOnEntities[MBTET][0].getOrder());

    int nb_gauss_pts;
    int rule = getRule(order);
    if(rule >= 0) {
      //if(mField.check_field(meshPositionsFieldName)) {
      //rule += 1;
      //}
      nb_gauss_pts = gm_rule_size(rule,3);
      gaussPts.resize(4,nb_gauss_pts);
      ierr = Grundmann_Moeller_integration_points_3D_TET(
        rule,&gaussPts(0,0),&gaussPts(1,0),&gaussPts(2,0),&gaussPts(3,0)
      ); CHKERRQ(ierr);
    } else {
      ierr = setGaussPts(order); CHKERRQ(ierr);
      nb_gauss_pts = gaussPts.size2();
    }

    ierr = shapeTETFunctions_H1(dataH1,&gaussPts(0,0),&gaussPts(1,0),&gaussPts(2,0),nb_gauss_pts); CHKERRQ(ierr);

    if((dataH1.spacesOnEntities[MBTRI]).test(HDIV)) {
      ierr = shapeTETFunctions_Hdiv(dataHdiv,&gaussPts(0,0),&gaussPts(1,0),&gaussPts(2,0),nb_gauss_pts); CHKERRQ(ierr);
    }

    if((dataH1.spacesOnEntities[MBTET]).test(L2)) {
      ierr = shapeTETFunctions_L2(dataL2,&gaussPts(0,0),&gaussPts(1,0),&gaussPts(2,0),nb_gauss_pts); CHKERRQ(ierr);
    }

    EntityHandle ent = fePtr->get_ent();
    int num_nodes;
    const EntityHandle* conn;
    rval = mField.get_moab().get_connectivity(ent,conn,num_nodes,true); CHKERR_PETSC(rval);
    coords.resize(num_nodes*3);
    rval = mField.get_moab().get_coords(conn,num_nodes,&*coords.data().begin()); CHKERR_PETSC(rval);
    vOlume = ShapeVolumeMBTET(&*dataH1.dataOnEntities[MBVERTEX][0].getDiffN().data().begin(),&*coords.data().begin());
    Jac.resize(3,3);
    invJac.resize(3,3);
    ierr = ShapeJacMBTET(
      &*dataH1.dataOnEntities[MBVERTEX][0].getDiffN().data().begin(),&*coords.begin(),&*Jac.data().begin()
    ); CHKERRQ(ierr);
    noalias(invJac) = Jac;
    ierr = ShapeInvJacMBTET(&*invJac.data().begin()); CHKERRQ(ierr);

    coordsAtGaussPts.resize(nb_gauss_pts,3);
    for(int gg = 0;gg<nb_gauss_pts;gg++) {
      for(int dd = 0;dd<3;dd++) {
        coordsAtGaussPts(gg,dd) = cblas_ddot(4,&dataH1.dataOnEntities[MBVERTEX][0].getN()(gg,0),1,&coords[dd],3);
      }
    }

    try {
      ierr = opSetInvJacH1.opRhs(dataH1); CHKERRQ(ierr);
      if((dataH1.spacesOnEntities[MBTET]).test(L2)) {
        ierr = opSetInvJacH1.opRhs(dataL2); CHKERRQ(ierr);
      }
      if((dataH1.spacesOnEntities[MBTRI]).test(HDIV)) {
        ierr = opPiolaTransform.opRhs(dataHdiv); CHKERRQ(ierr);
        ierr = opSetInvJacHdiv.opRhs(dataHdiv); CHKERRQ(ierr);
      }
    } catch (exception& ex) {
      ostringstream ss;
      ss << "thorw in method: " << ex.what() << " at line " << __LINE__ << " in file " << __FILE__;
      SETERRQ(PETSC_COMM_SELF,MOFEM_STD_EXCEPTION_THROW,ss.str().c_str());
    }

    if(mField.check_field(meshPositionsFieldName)) {
      BitFieldId id = mField.get_field_structure(meshPositionsFieldName)->get_id();
      if((fePtr->get_BitFieldId_data()&id).none()) {
        SETERRQ(PETSC_COMM_SELF,MOFEM_NOT_FOUND,"no MESH_NODE_POSITIONS in element data");
      }
      ierr = getEdgesOrder(dataH1,meshPositionsFieldName); CHKERRQ(ierr);
      ierr = getTrisOrder(dataH1,meshPositionsFieldName); CHKERRQ(ierr);
      ierr = getTetsOrder(dataH1,meshPositionsFieldName); CHKERRQ(ierr);
      ierr = getNodesFieldData(dataH1,meshPositionsFieldName); CHKERRQ(ierr);
      if(dataH1.dataOnEntities[MBVERTEX][0].getFieldData().size()!=12) {
        SETERRQ(PETSC_COMM_SELF,MOFEM_NOT_FOUND,"no MESH_NODE_POSITIONS in element data");
      }
      ierr = getEdgesFieldData(dataH1,meshPositionsFieldName); CHKERRQ(ierr);
      ierr = getTrisFieldData(dataH1,meshPositionsFieldName); CHKERRQ(ierr);
      ierr = getTetsFieldData(dataH1,meshPositionsFieldName); CHKERRQ(ierr);
      ierr = getNodesFieldDofs(dataH1,meshPositionsFieldName); CHKERRQ(ierr);
      ierr = getEdgesFieldDofs(dataH1,meshPositionsFieldName); CHKERRQ(ierr);
      ierr = getTrisFieldDofs(dataH1,meshPositionsFieldName); CHKERRQ(ierr);
      ierr = getTetsFieldDofs(dataH1,meshPositionsFieldName); CHKERRQ(ierr);
      try {
        ierr = opHOatGaussPoints.opRhs(dataH1); CHKERRQ(ierr);
        hoGaussPtsInvJac.resize(hoGaussPtsJac.size1(),hoGaussPtsJac.size2());
        ublas::noalias(hoGaussPtsInvJac) = hoGaussPtsJac;
        ublas::matrix<double> jac(3,3);
        hoGaussPtsDetJac.resize(nb_gauss_pts);
        for(int gg = 0;gg<nb_gauss_pts;gg++) {
          cblas_dcopy(9,&hoGaussPtsJac(gg,0),1,&jac(0,0),1);
          hoGaussPtsDetJac[gg] = ShapeDetJacMBTET(&jac(0,0));
          ierr = ShapeInvJacMBTET(&hoGaussPtsInvJac(gg,0)); CHKERRQ(ierr);
        }
        ierr = opSetHoInvJacH1.opRhs(dataH1); CHKERRQ(ierr);
        if((dataH1.spacesOnEntities[MBTET]).test(L2)) {
          ierr = opSetHoInvJacH1.opRhs(dataL2); CHKERRQ(ierr);
        }
        if((dataH1.spacesOnEntities[MBTRI]).test(HDIV)) {
          ierr = opSetHoPiolaTransform.opRhs(dataHdiv); CHKERRQ(ierr);
          ierr = opSetHoInvJacHdiv.opRhs(dataHdiv); CHKERRQ(ierr);
        }
      } catch (exception& ex) {
        ostringstream ss;
        ss << "problem with indices in method: " << ex.what() << " at line " << __LINE__ << " in file " << __FILE__;
        SETERRQ(PETSC_COMM_SELF,MOFEM_STD_EXCEPTION_THROW,ss.str().c_str());
      }
    } else {
      ublas::matrix<double> diffN(nb_gauss_pts,12);
      for(int gg = 0;gg<nb_gauss_pts;gg++) {
        for(int nn = 0;nn<4;nn++) {
          for(int dd = 0;dd<3;dd++) {
            diffN(gg,nn*3+dd) = dataH1.dataOnEntities[MBVERTEX][0].getDiffN()(nn,dd);
          }
        }
      }
      dataH1.dataOnEntities[MBVERTEX][0].getDiffN().resize(diffN.size1(),diffN.size2());
      dataH1.dataOnEntities[MBVERTEX][0].getDiffN().data().swap(diffN.data());
    }

    const UserDataOperator::OpType types[2] = {
      UserDataOperator::OPROW, UserDataOperator::OPCOL
    };
    vector<string> last_eval_field_name(2);
    DataForcesAndSurcesCore *op_data[2];
    FieldSpace space[2];

    boost::ptr_vector<UserDataOperator>::iterator oit,hi_oit;
    oit = opPtrVector.begin();
    hi_oit = opPtrVector.end();

    for(;oit!=hi_oit;oit++) {

      oit->setPtrFE(this);

      for(int ss = 0;ss!=2;ss++) {

        string field_name = !ss ? oit->rowFieldName : oit->colFieldName;
        BitFieldId data_id = mField.get_field_structure(field_name)->get_id();
        if((oit->getMoFEMFEPtr()->get_BitFieldId_data()&data_id).none()) {
          SETERRQ2(
            PETSC_COMM_SELF,MOFEM_DATA_INCONSISTENCT,"no data field < %s > on finite element < %s >",
            field_name.c_str(),feName.c_str()
          );
        }

        if(oit->getOpType()&types[ss] || oit->getOpType()&UserDataOperator::OPROWCOL) {

          space[ss] = mField.get_field_structure(field_name)->get_space();

          switch(space[ss]) {
            case H1:
            op_data[ss] = !ss ? &dataH1 : &derivedDataH1;
            break;
            case HCURL:
            SETERRQ(PETSC_COMM_SELF,MOFEM_NOT_IMPLEMENTED,"not implemented yet");
            break;
            case HDIV:
            op_data[ss] = !ss ? &dataHdiv : &derivedDataHdiv;
            break;
            case L2:
            op_data[ss] = !ss ? &dataL2 : &derivedDataL2;
            break;
            case NOFIELD:
            op_data[ss] = !ss ? &dataNoField : &dataNoFieldCol;
            break;
            case LASTSPACE:
            SETERRQ(PETSC_COMM_SELF,MOFEM_DATA_INCONSISTENCT,"unknown space");
            break;
          }

          if(last_eval_field_name[ss]!=field_name) {

            switch(space[ss]) {
              case H1:
              if(!ss) {
                ierr = getRowNodesIndices(*op_data[ss],field_name); CHKERRQ(ierr);
              } else {
                ierr = getColNodesIndices(*op_data[ss],field_name); CHKERRQ(ierr);
              }
              ierr = getNodesFieldData(*op_data[ss],field_name); CHKERRQ(ierr);
              ierr = getNodesFieldDofs(*op_data[ss],field_name); CHKERRQ(ierr);
              case HCURL:
              if(!ss) {
                ierr = getEdgesRowIndices(*op_data[ss],field_name); CHKERRQ(ierr);
              } else {
                ierr = getEdgesColIndices(*op_data[ss],field_name); CHKERRQ(ierr);
              }
              ierr = getEdgesOrder(*op_data[ss],field_name); CHKERRQ(ierr);
              ierr = getEdgesFieldData(*op_data[ss],field_name); CHKERRQ(ierr);
              ierr = getEdgesFieldDofs(*op_data[ss],field_name); CHKERRQ(ierr);
              case HDIV:
              if(!ss) {
                ierr = getTrisRowIndices(*op_data[ss],field_name); CHKERRQ(ierr);
              } else {
                ierr = getTrisColIndices(*op_data[ss],field_name); CHKERRQ(ierr);
              }
              ierr = getTrisOrder(*op_data[ss],field_name); CHKERRQ(ierr);
              ierr = getTrisFieldData(*op_data[ss],field_name); CHKERRQ(ierr);
              ierr = getTrisFieldDofs(*op_data[ss],field_name); CHKERRQ(ierr);
              case L2:
              if(!ss) {
                ierr = getTetsRowIndices(*op_data[ss],field_name); CHKERRQ(ierr);
              } else {

                ierr = getTetsColIndices(*op_data[ss],field_name); CHKERRQ(ierr);
              }
              ierr = getTetsOrder(*op_data[ss],field_name); CHKERRQ(ierr);
              ierr = getTetsFieldData(*op_data[ss],field_name); CHKERRQ(ierr);
              ierr = getTetsFieldDofs(*op_data[ss],field_name); CHKERRQ(ierr);
              break;
              case NOFIELD:
              if(!getNinTheLoop()) {
                // NOFIELD data arr the same for each element, can be retreived only once
                if(!ss) {
                  ierr = getNoFieldRowIndices(*op_data[ss],field_name); CHKERRQ(ierr);
                } else {
                  ierr = getNoFieldColIndices(*op_data[ss],field_name); CHKERRQ(ierr);
                }
                ierr = getNoFieldFieldData(*op_data[ss],field_name); CHKERRQ(ierr);
              }
              break;
              case LASTSPACE:
              SETERRQ(PETSC_COMM_SELF,MOFEM_DATA_INCONSISTENCT,"unknown space");
              break;
            }
            last_eval_field_name[ss]=field_name;

          }
        }
      }

      if(oit->getOpType()&UserDataOperator::OPROW) {
        try {
          ierr = oit->opRhs(*op_data[0]); CHKERRQ(ierr);
        } catch (exception& ex) {
          ostringstream ss;
          ss << "thorw in method: " << ex.what() << " at line " << __LINE__ << " in file " << __FILE__;
          SETERRQ(PETSC_COMM_SELF,MOFEM_STD_EXCEPTION_THROW,ss.str().c_str());
        }
      }


      if(oit->getOpType()&UserDataOperator::OPCOL) {
        try {
          ierr = oit->opRhs(*op_data[1]); CHKERRQ(ierr);
        } catch (exception& ex) {
          ostringstream ss;
          ss << "thorw in method: " << ex.what() << " at line " << __LINE__ << " in file " << __FILE__;
          SETERRQ(PETSC_COMM_SELF,MOFEM_STD_EXCEPTION_THROW,ss.str().c_str());
        }
      }


      if(oit->getOpType()&UserDataOperator::OPROWCOL) {
        try {
          ierr = oit->opLhs(*op_data[0],*op_data[1],oit->sYmm); CHKERRQ(ierr);
        } catch (exception& ex) {
          ostringstream ss;
          ss << "thorw in method: " << ex.what() << " at line " << __LINE__ << " in file " << __FILE__;
          SETERRQ(PETSC_COMM_SELF,MOFEM_STD_EXCEPTION_THROW,ss.str().c_str());
        }
      }

    }

  } catch (exception& ex) {
    ostringstream ss;
    ss << "thorw in method: " << ex.what() << " at line " << __LINE__ << " in file " << __FILE__;
    SETERRQ(PETSC_COMM_SELF,MOFEM_STD_EXCEPTION_THROW,ss.str().c_str());
  }

  PetscFunctionReturn(0);
}

PetscErrorCode VolumeElementForcesAndSourcesCore::UserDataOperator::getDivergenceMatrixOperato_Hdiv(
  int side,EntityType type,DataForcesAndSurcesCore::EntData &data,
  int gg,ublas::vector<FieldData> &div) {
    PetscFunctionBegin;

    try {

      int nb_dofs = data.getFieldData().size();
      if((unsigned int)nb_dofs != data.getDiffHdivN().size2()/9) {
        SETERRQ(PETSC_COMM_SELF,MOFEM_DATA_INCONSISTENCT,"data inconsistency");
      }

      if(nb_dofs == 0) PetscFunctionReturn(0);

      int dd = 0;
      for(;dd<nb_dofs;dd++) {
        div[dd] =
        (data.getDiffHdivN(dd,gg))(0,0)+
        (data.getDiffHdivN(dd,gg))(1,1)+
        (data.getDiffHdivN(dd,gg))(2,2);
      }

    } catch (exception& ex) {
      ostringstream ss;
      ss << "thorw in method: " << ex.what() << " at line " << __LINE__ << " in file " << __FILE__;
      SETERRQ(PETSC_COMM_SELF,1,ss.str().c_str());
    }

    PetscFunctionReturn(0);
  }

// **** Triangle ****

PetscErrorCode FaceElementForcesAndSourcesCore::operator()() {
  PetscFunctionBegin;

  if(fePtr->get_ent_type() != MBTRI) PetscFunctionReturn(0);

  ierr = getSpacesOnEntities(dataH1); CHKERRQ(ierr);

  //H1
  if((dataH1.spacesOnEntities[MBEDGE]).test(H1)) {
    ierr = getEdgesSense(dataH1); CHKERRQ(ierr);
    ierr = getTrisSense(dataH1); CHKERRQ(ierr);
    ierr = getEdgesOrder(dataH1,H1); CHKERRQ(ierr);
    ierr = getTrisOrder(dataH1,H1); CHKERRQ(ierr);
  }

  //Hdiv
  if((dataH1.spacesOnEntities[MBTRI]).test(HDIV)) {
    ierr = getTrisSense(dataHdiv); CHKERRQ(ierr);
    ierr = getTrisOrder(dataHdiv,HDIV); CHKERRQ(ierr);
  }

  int order = 1;
  for(unsigned int ee = 0;ee<dataH1.dataOnEntities[MBEDGE].size();ee++) {
    order = max(order,dataH1.dataOnEntities[MBEDGE][ee].getOrder());
  }
  for(unsigned int ff = 0;ff<dataHdiv.dataOnEntities[MBTRI].size();ff++) {
    order = max(order,dataHdiv.dataOnEntities[MBTRI][ff].getOrder());
  }

  int nb_gauss_pts;
  int rule = getRule(order);
  if(rule >= 0) {
    //if(mField.check_field(meshPositionsFieldName)) {
    //rule += 1;
    //}
    nb_gauss_pts = gm_rule_size(rule,2);
    gaussPts.resize(3,nb_gauss_pts);
    ierr = Grundmann_Moeller_integration_points_2D_TRI(
      rule,&gaussPts(0,0),&gaussPts(1,0),&gaussPts(2,0)
    ); CHKERRQ(ierr);
  } else {
    ierr = setGaussPts(order); CHKERRQ(ierr);
    nb_gauss_pts = gaussPts.size2();
  }

  ierr = shapeTRIFunctions_H1(dataH1,&gaussPts(0,0),&gaussPts(1,0),nb_gauss_pts); CHKERRQ(ierr);

  if(dataH1.spacesOnEntities[MBTRI].test(HDIV)) {
    ierr = shapeTRIFunctions_Hdiv(dataHdiv,&gaussPts(0,0),&gaussPts(1,0),nb_gauss_pts); CHKERRQ(ierr); CHKERRQ(ierr);
  }

  EntityHandle ent = fePtr->get_ent();
  int num_nodes;
  const EntityHandle* conn;
  rval = mField.get_moab().get_connectivity(ent,conn,num_nodes,true); CHKERR_PETSC(rval);
  coords.resize(num_nodes*3);
  rval = mField.get_moab().get_coords(conn,num_nodes,&*coords.data().begin()); CHKERR_PETSC(rval);

  normal.resize(3);
  ierr = ShapeFaceNormalMBTRI(
    &*dataH1.dataOnEntities[MBVERTEX][0].getDiffN().data().begin(),
    &*coords.data().begin(),&*normal.data().begin()
  ); CHKERRQ(ierr);
  aRea = cblas_dnrm2(3,&*normal.data().begin(),1)*0.5;

  coordsAtGaussPts.resize(nb_gauss_pts,3);
  for(int gg = 0;gg<nb_gauss_pts;gg++) {
    for(int dd = 0;dd<3;dd++) {
      coordsAtGaussPts(gg,dd) = cblas_ddot(3,&dataH1.dataOnEntities[MBVERTEX][0].getN()(gg,0),1,&coords[dd],3);
    }
  }

  // In linear geometry derivatives are constant,
  // this in expense of efficiency makes implementation
  // constant between vertices and other types of entities
  ublas::matrix<double> diffN(nb_gauss_pts,6);
  for(int gg = 0;gg<nb_gauss_pts;gg++) {
    for(int nn = 0;nn<3;nn++) {
      for(int dd = 0;dd<2;dd++) {
        diffN(gg,nn*2+dd) = dataH1.dataOnEntities[MBVERTEX][0].getDiffN()(nn,dd);
      }
    }
  }
  dataH1.dataOnEntities[MBVERTEX][0].getDiffN().resize(diffN.size1(),diffN.size2());
  dataH1.dataOnEntities[MBVERTEX][0].getDiffN().data().swap(diffN.data());

  if(mField.check_field(meshPositionsFieldName)) {
    nOrmals_at_GaussPt.resize(nb_gauss_pts,3);
    tAngent1_at_GaussPt.resize(nb_gauss_pts,3);
    tAngent2_at_GaussPt.resize(nb_gauss_pts,3);
    ierr = getEdgesOrder(dataH1,meshPositionsFieldName); CHKERRQ(ierr);
    ierr = getTrisOrder(dataH1,meshPositionsFieldName); CHKERRQ(ierr);
    ierr = getNodesFieldData(dataH1,meshPositionsFieldName); CHKERRQ(ierr);
    ierr = getEdgesFieldData(dataH1,meshPositionsFieldName); CHKERRQ(ierr);
    ierr = getTrisFieldData(dataH1,meshPositionsFieldName); CHKERRQ(ierr);
    ierr = getNodesFieldDofs(dataH1,meshPositionsFieldName); CHKERRQ(ierr);
    ierr = getEdgesFieldDofs(dataH1,meshPositionsFieldName); CHKERRQ(ierr);
    ierr = getTrisFieldDofs(dataH1,meshPositionsFieldName); CHKERRQ(ierr);
    try {
      ierr = opHONormals.opRhs(dataH1); CHKERRQ(ierr);
      ierr = opHONormals.calculateNormals(); CHKERRQ(ierr);
    } catch (exception& ex) {
      ostringstream ss;
      ss << "thorw in method: " << ex.what() << " at line " << __LINE__ << " in file " << __FILE__;
      SETERRQ(PETSC_COMM_SELF,1,ss.str().c_str());
    }
  }

  if(dataH1.spacesOnEntities[MBTRI].test(HDIV)) {
    ierr = opSetPiolaTransoformOnTriangle.opRhs(dataHdiv); CHKERRQ(ierr);
  }

  const UserDataOperator::OpType types[2] = {
    UserDataOperator::OPROW, UserDataOperator::OPCOL
  };
  vector<string> last_eval_field_name(2);
  DataForcesAndSurcesCore *op_data[2];
  FieldSpace space[2];

  boost::ptr_vector<UserDataOperator>::iterator oit,hi_oit;
  oit = opPtrVector.begin();
  hi_oit = opPtrVector.end();

  for(;oit!=hi_oit;oit++) {

    oit->setPtrFE(this);

    for(int ss = 0;ss!=2;ss++) {

      string field_name = !ss ? oit->rowFieldName : oit->colFieldName;
      BitFieldId data_id = mField.get_field_structure(field_name)->get_id();
      if((oit->getMoFEMFEPtr()->get_BitFieldId_data()&data_id).none()) {
        SETERRQ2(
          PETSC_COMM_SELF,MOFEM_DATA_INCONSISTENCT,"no data field < %s > on finite element < %s >",
          field_name.c_str(),feName.c_str()
        );
      }

      if(oit->getOpType()&types[ss] || oit->getOpType()&UserDataOperator::OPROWCOL) {

        space[ss] = mField.get_field_structure(field_name)->get_space();

        switch(space[ss]) {
          case H1:
          op_data[ss] = !ss ? &dataH1 : &derivedDataH1;
          break;
          case HCURL:
          SETERRQ(PETSC_COMM_SELF,MOFEM_NOT_IMPLEMENTED,"not implemented yet");
          break;
          case HDIV:
          op_data[ss] = !ss ? &dataHdiv : &derivedDataHdiv;
          break;
          case L2:
          SETERRQ(PETSC_COMM_SELF,MOFEM_NOT_IMPLEMENTED,"not make sanes on face");
          break;
          case NOFIELD:
          op_data[ss] = !ss ? &dataNoField : &dataNoFieldCol;
          break;
          case LASTSPACE:
          SETERRQ(PETSC_COMM_SELF,MOFEM_DATA_INCONSISTENCT,"unknown space");
          break;
        }

        if(last_eval_field_name[ss]!=field_name) {

          switch(space[ss]) {
            case H1:
            if(!ss) {
              ierr = getRowNodesIndices(*op_data[ss],field_name); CHKERRQ(ierr);
            } else {
              ierr = getColNodesIndices(*op_data[ss],field_name); CHKERRQ(ierr);
            }
            ierr = getNodesFieldData(*op_data[ss],field_name); CHKERRQ(ierr);
            ierr = getNodesFieldDofs(*op_data[ss],field_name); CHKERRQ(ierr);
            case HCURL:
            if(!ss) {
              ierr = getEdgesRowIndices(*op_data[ss],field_name); CHKERRQ(ierr);
            } else {
              ierr = getEdgesColIndices(*op_data[ss],field_name); CHKERRQ(ierr);
            }
            ierr = getEdgesOrder(*op_data[ss],field_name); CHKERRQ(ierr);
            ierr = getEdgesFieldData(*op_data[ss],field_name); CHKERRQ(ierr);
            ierr = getEdgesFieldDofs(*op_data[ss],field_name); CHKERRQ(ierr);
            case HDIV:
            if(!ss) {
              ierr = getTrisRowIndices(*op_data[ss],field_name); CHKERRQ(ierr);
            } else {
              ierr = getTrisColIndices(*op_data[ss],field_name); CHKERRQ(ierr);
            }
            ierr = getTrisOrder(*op_data[ss],field_name); CHKERRQ(ierr);
            ierr = getTrisFieldData(*op_data[ss],field_name); CHKERRQ(ierr);
            ierr = getTrisFieldDofs(*op_data[ss],field_name); CHKERRQ(ierr);
            break;
            case L2:
            SETERRQ(PETSC_COMM_SELF,MOFEM_NOT_IMPLEMENTED,"not make sanes on face");
            break;
            case NOFIELD:
            if(!getNinTheLoop()) {
              // NOFIELD data arr the same for each element, can be retreived only once
              if(!ss) {
                ierr = getNoFieldRowIndices(*op_data[ss],field_name); CHKERRQ(ierr);
              } else {
                ierr = getNoFieldColIndices(*op_data[ss],field_name); CHKERRQ(ierr);
              }
              ierr = getNoFieldFieldData(*op_data[ss],field_name); CHKERRQ(ierr);
            }
            break;
            case LASTSPACE:
            SETERRQ(PETSC_COMM_SELF,MOFEM_DATA_INCONSISTENCT,"unknown space");
            break;
          }
          last_eval_field_name[ss]=field_name;

        }
      }
    }

    if(oit->getOpType()&UserDataOperator::OPROW) {
      try {
        ierr = oit->opRhs(*op_data[0]); CHKERRQ(ierr);
      } catch (exception& ex) {
        ostringstream ss;
        ss << "thorw in method: " << ex.what() << " at line " << __LINE__ << " in file " << __FILE__;
        SETERRQ(PETSC_COMM_SELF,MOFEM_STD_EXCEPTION_THROW,ss.str().c_str());
      }
    }


    if(oit->getOpType()&UserDataOperator::OPCOL) {
      try {
        ierr = oit->opRhs(*op_data[1]); CHKERRQ(ierr);
      } catch (exception& ex) {
        ostringstream ss;
        ss << "thorw in method: " << ex.what() << " at line " << __LINE__ << " in file " << __FILE__;
        SETERRQ(PETSC_COMM_SELF,MOFEM_STD_EXCEPTION_THROW,ss.str().c_str());
      }
    }


    if(oit->getOpType()&UserDataOperator::OPROWCOL) {
      try {
        ierr = oit->opLhs(*op_data[0],*op_data[1],oit->sYmm); CHKERRQ(ierr);
      } catch (exception& ex) {
        ostringstream ss;
        ss << "thorw in method: " << ex.what() << " at line " << __LINE__ << " in file " << __FILE__;
        SETERRQ(PETSC_COMM_SELF,MOFEM_STD_EXCEPTION_THROW,ss.str().c_str());
      }
    }

  }

  PetscFunctionReturn(0);
}

// **** Edge ****

PetscErrorCode EdgeElementForcesAndSurcesCore::operator()() {
  PetscFunctionBegin;

  if(fePtr->get_ent_type() != MBEDGE) PetscFunctionReturn(0);

  //PetscAttachDebugger();

  ierr = getEdgesOrder(dataH1,H1); CHKERRQ(ierr);

  int order = dataH1.dataOnEntities[MBEDGE][0].getOrder();
  int rule = getRule(order);
  int nb_gauss_pts = gm_rule_size(rule,1);
  gaussPts.resize(2,nb_gauss_pts);

  ierr = Grundmann_Moeller_integration_points_1D_EDGE(rule,&gaussPts(0,0),&gaussPts(1,0)); CHKERRQ(ierr);
  ierr = shapeEDGEFunctions_H1(dataH1,&gaussPts(0,0),nb_gauss_pts); CHKERRQ(ierr);

  EntityHandle ent = fePtr->get_ent();
  int num_nodes;
  const EntityHandle* conn;
  rval = mField.get_moab().get_connectivity(ent,conn,num_nodes,true); CHKERR_PETSC(rval);
  coords.resize(num_nodes*3);
  rval = mField.get_moab().get_coords(conn,num_nodes,&*coords.data().begin()); CHKERR_PETSC(rval);

  dIrection.resize(3);
  cblas_dcopy(3,&coords[3],1,&*dIrection.data().begin(),1);
  cblas_daxpy(3,-1.,&coords[0],1,&*dIrection.data().begin(),1);
  lEngth = cblas_dnrm2(3,&*dIrection.data().begin(),1);

  coordsAtGaussPts.resize(nb_gauss_pts,3);
  for(int gg = 0;gg<nb_gauss_pts;gg++) {
    for(int dd = 0;dd<3;dd++) {
      coordsAtGaussPts(gg,dd)
      = N_MBEDGE0(gaussPts(0,gg))*coords[dd] + N_MBEDGE1(gaussPts(0,gg))*coords[3+dd];
    }
  }
  //cerr << coordsAtGaussPts << endl;

  const UserDataOperator::OpType types[2] = {
    UserDataOperator::OPROW, UserDataOperator::OPCOL
  };
  vector<string> last_eval_field_name(2);
  DataForcesAndSurcesCore *op_data[2];
  FieldSpace space[2];

  boost::ptr_vector<UserDataOperator>::iterator oit,hi_oit;
  oit = opPtrVector.begin();
  hi_oit = opPtrVector.end();

  for(;oit!=hi_oit;oit++) {

    oit->setPtrFE(this);

    for(int ss = 0;ss!=2;ss++) {

      string field_name = !ss ? oit->rowFieldName : oit->colFieldName;
      BitFieldId data_id = mField.get_field_structure(field_name)->get_id();
      if((oit->getMoFEMFEPtr()->get_BitFieldId_data()&data_id).none()) {
        SETERRQ2(
          PETSC_COMM_SELF,MOFEM_DATA_INCONSISTENCT,"no data field < %s > on finite element < %s >",
          field_name.c_str(),feName.c_str()
        );
      }

      if(oit->getOpType()&types[ss] || oit->getOpType()&UserDataOperator::OPROWCOL) {

        space[ss] = mField.get_field_structure(field_name)->get_space();

        switch(space[ss]) {
          case H1:
          op_data[ss] = !ss ? &dataH1 : &derivedDataH1;
          break;
          case HCURL:
          SETERRQ(PETSC_COMM_SELF,MOFEM_NOT_IMPLEMENTED,"not implemented yet");
          break;
          case HDIV:
          SETERRQ(PETSC_COMM_SELF,MOFEM_NOT_IMPLEMENTED,"not make sanes on edge");
          break;
          case L2:
          SETERRQ(PETSC_COMM_SELF,MOFEM_NOT_IMPLEMENTED,"not make sanes on edge");
          break;
          case NOFIELD:
          op_data[ss] = !ss ? &dataNoField : &dataNoFieldCol;
          break;
          case LASTSPACE:
          SETERRQ(PETSC_COMM_SELF,MOFEM_DATA_INCONSISTENCT,"unknown space");
          break;
        }

        if(last_eval_field_name[ss]!=field_name) {

          switch(space[ss]) {
            case H1:
            if(!ss) {
              ierr = getRowNodesIndices(*op_data[ss],field_name); CHKERRQ(ierr);
            } else {
              ierr = getColNodesIndices(*op_data[ss],field_name); CHKERRQ(ierr);
            }
            ierr = getNodesFieldData(*op_data[ss],field_name); CHKERRQ(ierr);
            ierr = getNodesFieldDofs(*op_data[ss],field_name); CHKERRQ(ierr);
            case HCURL:
            if(!ss) {
              ierr = getEdgesRowIndices(*op_data[ss],field_name); CHKERRQ(ierr);
            } else {
              ierr = getEdgesColIndices(*op_data[ss],field_name); CHKERRQ(ierr);
            }
            ierr = getEdgesOrder(*op_data[ss],field_name); CHKERRQ(ierr);
            ierr = getEdgesFieldData(*op_data[ss],field_name); CHKERRQ(ierr);
            ierr = getEdgesFieldDofs(*op_data[ss],field_name); CHKERRQ(ierr);
            break;
            case HDIV:
            SETERRQ(PETSC_COMM_SELF,MOFEM_NOT_IMPLEMENTED,"not make sanes on edge");
            break;
            case L2:
            SETERRQ(PETSC_COMM_SELF,MOFEM_NOT_IMPLEMENTED,"not make sanes on edge");
            break;
            case NOFIELD:
            if(!getNinTheLoop()) {
              // NOFIELD data arr the same for each element, can be retreived only once
              if(!ss) {
                ierr = getNoFieldRowIndices(*op_data[ss],field_name); CHKERRQ(ierr);
              } else {
                ierr = getNoFieldColIndices(*op_data[ss],field_name); CHKERRQ(ierr);
              }
              ierr = getNoFieldFieldData(*op_data[ss],field_name); CHKERRQ(ierr);
            }
            break;
            case LASTSPACE:
            SETERRQ(PETSC_COMM_SELF,MOFEM_DATA_INCONSISTENCT,"unknown space");
            break;
          }
          last_eval_field_name[ss]=field_name;

        }
      }
    }

    if(oit->getOpType()&UserDataOperator::OPROW) {
      try {
        ierr = oit->opRhs(*op_data[0]); CHKERRQ(ierr);
      } catch (exception& ex) {
        ostringstream ss;
        ss << "thorw in method: " << ex.what() << " at line " << __LINE__ << " in file " << __FILE__;
        SETERRQ(PETSC_COMM_SELF,MOFEM_STD_EXCEPTION_THROW,ss.str().c_str());
      }
    }


    if(oit->getOpType()&UserDataOperator::OPCOL) {
      try {
        ierr = oit->opRhs(*op_data[1]); CHKERRQ(ierr);
      } catch (exception& ex) {
        ostringstream ss;
        ss << "thorw in method: " << ex.what() << " at line " << __LINE__ << " in file " << __FILE__;
        SETERRQ(PETSC_COMM_SELF,MOFEM_STD_EXCEPTION_THROW,ss.str().c_str());
      }
    }


    if(oit->getOpType()&UserDataOperator::OPROWCOL) {
      try {
        ierr = oit->opLhs(*op_data[0],*op_data[1],oit->sYmm); CHKERRQ(ierr);
      } catch (exception& ex) {
        ostringstream ss;
        ss << "thorw in method: " << ex.what() << " at line " << __LINE__ << " in file " << __FILE__;
        SETERRQ(PETSC_COMM_SELF,MOFEM_STD_EXCEPTION_THROW,ss.str().c_str());
      }
    }

  }


  /*string last_row_eval_field_name;
  string last_col_eval_field_name;

  for(int ss = 0;ss<2;ss++) {

    boost::ptr_vector<ForcesAndSurcesCore::UserDataOperator>::iterator oit,hi_oit;
    if(!ss) {

      oit = rowOpPtrVector.begin();
      hi_oit = rowOpPtrVector.end();

    } else {

      oit = colOpPtrVector.begin();
      hi_oit = colOpPtrVector.end();

    }

    for(;oit!=hi_oit;oit++) {

      oit->setPtrFE(this);

      BitFieldId data_id = mField.get_field_structure(oit->rowFieldName)->get_id();
      if((oit->getMoFEMFEPtr()->get_BitFieldId_data()&data_id).none()) {
        SETERRQ2(
          PETSC_COMM_SELF,MOFEM_DATA_INCONSISTENCT,"no data field < %s > on finite element < %s >",
          oit->rowFieldName.c_str(),feName.c_str()
        );
      }

      FieldSpace space;
      if(!ss) {
        space = mField.get_field_structure(oit->rowFieldName)->get_space();
      } else {
        space = mField.get_field_structure(oit->colFieldName)->get_space();
      }

      DataForcesAndSurcesCore *op_data = NULL;
      switch(space) {
        case H1:
        op_data = &dataH1;
        break;
        case HCURL:
        SETERRQ(PETSC_COMM_SELF,MOFEM_NOT_IMPLEMENTED,"not implemented yet");
        break;
        case HDIV:
        break;
        case L2:
        break;
        default:
        SETERRQ(PETSC_COMM_SELF,MOFEM_DATA_INCONSISTENCT,"data inconsistency");
        break;
      }

      string data_field_name;
      bool get_data = false;

      if(!ss) {

        data_field_name = oit->rowFieldName;

        if(last_row_eval_field_name!=oit->rowFieldName) {

          get_data = true;

          switch(space) {
            case H1:
            ierr = getRowNodesIndices(*op_data,oit->rowFieldName); CHKERRQ(ierr);
            case HCURL:
            ierr = getEdgesRowIndices(*op_data,oit->rowFieldName); CHKERRQ(ierr);
            default:
            break;
          }

          last_row_eval_field_name=oit->rowFieldName;

        }

      } else {

        data_field_name = oit->rowFieldName;


        if(last_col_eval_field_name!=oit->colFieldName) {

          get_data = true;

          switch(space) {
            case H1:
            ierr = getColNodesIndices(*op_data,oit->colFieldName); CHKERRQ(ierr);
            case HCURL:
            ierr = getEdgesColIndices(*op_data,oit->colFieldName); CHKERRQ(ierr);
            default:
            break;
          }

          last_col_eval_field_name=oit->colFieldName;

        }

      }

      if(get_data) {

        switch(space) {
          case H1:
          ierr = getNodesFieldData(*op_data,data_field_name); CHKERRQ(ierr);
          ierr = getNodesFieldDofs(*op_data,data_field_name); CHKERRQ(ierr);
          case HCURL:
          ierr = getEdgesOrder(*op_data,data_field_name); CHKERRQ(ierr);
          ierr = getEdgesFieldData(*op_data,data_field_name); CHKERRQ(ierr);
          ierr = getEdgesFieldDofs(*op_data,data_field_name); CHKERRQ(ierr);
          default:
          break;
        }

      }

      try {
        ierr = oit->opRhs(*op_data); CHKERRQ(ierr);
      } catch (exception& ex) {
        ostringstream ss;
        ss << "thorw in method: " << ex.what() << " at line " << __LINE__ << " in file " << __FILE__;
        SETERRQ(PETSC_COMM_SELF,MOFEM_STD_EXCEPTION_THROW,ss.str().c_str());
      }

    }
  }

  last_row_eval_field_name.clear();
  last_col_eval_field_name.clear();

  for(
    boost::ptr_vector<ForcesAndSurcesCore::UserDataOperator>::iterator oit = rowColOpPtrVector.begin();
    oit != rowColOpPtrVector.end(); oit++
  ) {

    oit->setPtrFE(this);

    BitFieldId row_id = mField.get_field_structure(oit->rowFieldName)->get_id();
    BitFieldId col_id = mField.get_field_structure(oit->colFieldName)->get_id();

    if((oit->getMoFEMFEPtr()->get_BitFieldId_data()&row_id).none()) {
      SETERRQ1(PETSC_COMM_SELF,MOFEM_NOT_FOUND,"no row field < %s > on finite elemeny",oit->rowFieldName.c_str());
    }
    if((oit->getMoFEMFEPtr()->get_BitFieldId_data()&col_id).none()) {
      SETERRQ1(PETSC_COMM_SELF,MOFEM_NOT_FOUND,"no data field < %s > on finite elemeny",oit->colFieldName.c_str());
    }

    FieldSpace row_space = mField.get_field_structure(oit->rowFieldName)->get_space();

    DataForcesAndSurcesCore *row_op_data = NULL;
    switch(row_space) {
      case H1:
      row_op_data = &dataH1;
      break;
      case HCURL:
      SETERRQ(PETSC_COMM_SELF,MOFEM_NOT_IMPLEMENTED,"not implemented yet");
      break;
      case HDIV:
      case L2:
      default:
      SETERRQ(PETSC_COMM_SELF,MOFEM_DATA_INCONSISTENCT,"data inconsistency");
      break;
    }

    if(last_row_eval_field_name!=oit->rowFieldName) {

      switch(row_space) {
        case H1:
        ierr = getRowNodesIndices(*row_op_data,oit->rowFieldName); CHKERRQ(ierr);
        ierr = getNodesFieldData(*row_op_data,oit->rowFieldName); CHKERRQ(ierr);
        ierr = getNodesFieldDofs(*row_op_data,oit->rowFieldName); CHKERRQ(ierr);
        case HCURL:
        ierr = getEdgesRowIndices(*row_op_data,oit->rowFieldName); CHKERRQ(ierr);
        ierr = getEdgesOrder(*row_op_data,oit->rowFieldName); CHKERRQ(ierr);
        ierr = getEdgesFieldData(*row_op_data,oit->rowFieldName); CHKERRQ(ierr);
        ierr = getEdgesFieldDofs(*row_op_data,oit->rowFieldName); CHKERRQ(ierr);
        default:
        break;
      }

      last_row_eval_field_name=oit->rowFieldName;

    }

    FieldSpace col_space = mField.get_field_structure(oit->colFieldName)->get_space();
    DataForcesAndSurcesCore *col_op_data = NULL;
    switch(col_space) {
      case H1:
      col_op_data = &derivedDataH1;
      break;
      case HCURL:
      SETERRQ(PETSC_COMM_SELF,MOFEM_NOT_IMPLEMENTED,"not implemented yet");
      break;
      case HDIV:
      case L2:
      default:
      SETERRQ(PETSC_COMM_SELF,MOFEM_DATA_INCONSISTENCT,"data inconsistency");
      break;
    }

    if(last_col_eval_field_name!=oit->colFieldName) {

      switch(col_space) {
        case H1:
        ierr = getColNodesIndices(*col_op_data,oit->colFieldName); CHKERRQ(ierr);
        ierr = getNodesFieldData(*col_op_data,oit->colFieldName); CHKERRQ(ierr);
        ierr = getNodesFieldDofs(*col_op_data,oit->colFieldName); CHKERRQ(ierr);
        case HCURL:
        ierr = getEdgesColIndices(*col_op_data,oit->colFieldName); CHKERRQ(ierr);
        ierr = getEdgesOrder(*col_op_data,oit->colFieldName); CHKERRQ(ierr);
        ierr = getEdgesFieldData(*col_op_data,oit->colFieldName); CHKERRQ(ierr);
        ierr = getEdgesFieldDofs(*col_op_data,oit->colFieldName); CHKERRQ(ierr);
        default:
        break;
      }

      last_col_eval_field_name=oit->colFieldName;

    }

    try {
      ierr = oit->opLhs(*row_op_data,*col_op_data,oit->sYmm); CHKERRQ(ierr);
    } catch (exception& ex) {
      ostringstream ss;
      ss << "thorw in method: " << ex.what() << " at line " << __LINE__ << " in file " << __FILE__;
      SETERRQ(PETSC_COMM_SELF,MOFEM_STD_EXCEPTION_THROW,ss.str().c_str());
    }

  }*/

  PetscFunctionReturn(0);
}

// Vertex

PetscErrorCode VertexElementForcesAndSourcesCore::operator()() {
  PetscFunctionBegin;

  if(fePtr->get_ent_type() != MBVERTEX) PetscFunctionReturn(0);

  EntityHandle ent = fePtr->get_ent();
  coords.resize(3);
  rval = mField.get_moab().get_coords(&ent,1,&*coords.data().begin()); CHKERR_PETSC(rval);

  /*DataForcesAndSurcesCore *col_data = &derivedData;

  for(
    boost::ptr_vector<ForcesAndSurcesCore::UserDataOperator>::iterator oit = rowOpPtrVector.begin();
    oit != rowOpPtrVector.end(); oit++
  ) {

    oit->setPtrFE(this);
    BitFieldId row_id = mField.get_field_structure(oit->rowFieldName)->get_id();
    BitFieldId col_id = mField.get_field_structure(oit->colFieldName)->get_id();

    if((oit->getMoFEMFEPtr()->get_BitFieldId_row()&row_id).none()) {
      SETERRQ1(PETSC_COMM_SELF,1,"no row field < %s > on finite elemeny",oit->rowFieldName.c_str());
    }
    if((oit->getMoFEMFEPtr()->get_BitFieldId_data()&col_id).none()) {
      SETERRQ1(PETSC_COMM_SELF,1,"no data field < %s > on finite elemeny",oit->rowFieldName.c_str());
    }

    ierr = getRowNodesIndices(data,oit->rowFieldName); CHKERRQ(ierr);
    ierr = getNodesFieldData(data,oit->colFieldName); CHKERRQ(ierr);
    ierr = getNodesFieldDofs(data,oit->colFieldName); CHKERRQ(ierr);

    try {
      ierr = oit->opRhs(data); CHKERRQ(ierr);
    } catch (exception& ex) {
      ostringstream ss;
      ss << "thorw in method: " << ex.what() << " at line " << __LINE__ << " in file " << __FILE__;
      SETERRQ(PETSC_COMM_SELF,1,ss.str().c_str());
    }

  }

  for(
    boost::ptr_vector<ForcesAndSurcesCore::UserDataOperator>::iterator oit = rowColOpPtrVector.begin();
    oit != rowColOpPtrVector.end(); oit++
  ) {

    oit->setPtrFE(this);
    BitFieldId row_id = mField.get_field_structure(oit->rowFieldName)->get_id();
    BitFieldId col_id = mField.get_field_structure(oit->colFieldName)->get_id();

    if((oit->getMoFEMFEPtr()->get_BitFieldId_row()&row_id).none()) {
      SETERRQ1(PETSC_COMM_SELF,1,"no row field < %s > on finite elemeny",oit->rowFieldName.c_str());
    }
    if((oit->getMoFEMFEPtr()->get_BitFieldId_col()&col_id).none()) {
      SETERRQ1(PETSC_COMM_SELF,1,"no data field < %s > on finite elemeny",oit->colFieldName.c_str());
    }
    if((oit->getMoFEMFEPtr()->get_BitFieldId_data()&col_id).none()) {
      SETERRQ1(PETSC_COMM_SELF,1,"no data field < %s > on finite elemeny",oit->colFieldName.c_str());
    }

    ierr = getRowNodesIndices(data,oit->rowFieldName); CHKERRQ(ierr);
    ierr = getColNodesIndices(*col_data,oit->colFieldName); CHKERRQ(ierr);
    ierr = getNodesFieldData(*col_data,oit->colFieldName); CHKERRQ(ierr);
    ierr = getNodesFieldDofs(*col_data,oit->colFieldName); CHKERRQ(ierr);

    try {
      ierr = oit->opLhs(data,*col_data,true); CHKERRQ(ierr);
    } catch (exception& ex) {
      ostringstream ss;
      ss << "thorw in method: " << ex.what() << " at line " << __LINE__ << " in file " << __FILE__;
      SETERRQ(PETSC_COMM_SELF,1,ss.str().c_str());
    }

  }*/

  PetscFunctionReturn(0);
}

PetscErrorCode FlatPrismElementForcesAndSurcesCore::operator()() {
  PetscFunctionBegin;

  if(fePtr->get_ent_type() != MBPRISM) PetscFunctionReturn(0);

  try {

    ierr = getSpacesOnEntities(dataH1); CHKERRQ(ierr);

    //H1
    if((dataH1.spacesOnEntities[MBEDGE]).test(H1)) {
      ierr = getEdgesSense(dataH1); CHKERRQ(ierr);
      ierr = getTrisSense(dataH1); CHKERRQ(ierr);
      ierr = getEdgesOrder(dataH1,H1); CHKERRQ(ierr);
      ierr = getTrisOrder(dataH1,H1); CHKERRQ(ierr);
    }

    //Hdiv
    if((dataH1.spacesOnEntities[MBTRI]).test(HDIV)) {
      ierr = getTrisSense(dataHdiv); CHKERRQ(ierr);
      ierr = getTrisOrder(dataHdiv,HDIV); CHKERRQ(ierr);
    }

    int order = 1;
    for(unsigned int ee = 0;ee<dataH1.dataOnEntities[MBEDGE].size();ee++) {
      order = max(order,dataH1.dataOnEntities[MBEDGE][ee].getOrder());
    }
    for(unsigned int ff = 0;ff<dataHdiv.dataOnEntities[MBTRI].size();ff++) {
      order = max(order,dataHdiv.dataOnEntities[MBTRI][ff].getOrder());
    }

    int nb_gauss_pts;
    int rule = getRule(order);
    if(rule >= 0) {
      nb_gauss_pts = gm_rule_size(rule,2);
      gaussPts.resize(3,nb_gauss_pts);
      ierr = Grundmann_Moeller_integration_points_2D_TRI(
        rule,&gaussPts(0,0),&gaussPts(1,0),&gaussPts(2,0)
      ); CHKERRQ(ierr);
    } else {
      ierr = setGaussPts(order); CHKERRQ(ierr);
      nb_gauss_pts = gaussPts.size2();
    }

    ierr = shapeFlatPRISMFunctions_H1(dataH1,&gaussPts(0,0),&gaussPts(1,0),nb_gauss_pts); CHKERRQ(ierr);
    if(dataH1.spacesOnEntities[MBTRI].test(HDIV)) {
      ierr = shapeFlatPRISMFunctions_Hdiv(dataHdiv,&gaussPts(0,0),&gaussPts(1,0),nb_gauss_pts); CHKERRQ(ierr); CHKERRQ(ierr);
    }

    EntityHandle ent = fePtr->get_ent();
    int num_nodes;
    const EntityHandle* conn;
    rval = mField.get_moab().get_connectivity(ent,conn,num_nodes,true); CHKERR_PETSC(rval);
    coords.resize(num_nodes*3);
    rval = mField.get_moab().get_coords(conn,num_nodes,&*coords.data().begin()); CHKERR_PETSC(rval);

    normal.resize(3);
    ierr = ShapeFaceNormalMBTRI(
      &*dataH1.dataOnEntities[MBVERTEX][0].getDiffN().data().begin(),
      &*coords.data().begin(),&*normal.data().begin()
    ); CHKERRQ(ierr);
    aRea = cblas_dnrm2(3,&*normal.data().begin(),1)*0.5;

    coordsAtGaussPts.resize(nb_gauss_pts,3);
    for(int gg = 0;gg<nb_gauss_pts;gg++) {
      for(int dd = 0;dd<3;dd++) {
        coordsAtGaussPts(gg,dd) = cblas_ddot(3,&dataH1.dataOnEntities[MBVERTEX][0].getN()(gg,0),1,&coords[dd],3);
      }
    }

    try {
      if(mField.check_field(meshPositionsFieldName)) {
        nOrmals_at_GaussPtF3.resize(nb_gauss_pts,3);
        tAngent1_at_GaussPtF3.resize(nb_gauss_pts,3);
        tAngent2_at_GaussPtF3.resize(nb_gauss_pts,3);
        nOrmals_at_GaussPtF4.resize(nb_gauss_pts,3);
        tAngent1_at_GaussPtF4.resize(nb_gauss_pts,3);
        tAngent2_at_GaussPtF4.resize(nb_gauss_pts,3);
        ierr = getEdgesOrder(dataH1,meshPositionsFieldName); CHKERRQ(ierr);
        ierr = getTrisOrder(dataH1,meshPositionsFieldName); CHKERRQ(ierr);
        ierr = getNodesFieldData(dataH1,meshPositionsFieldName); CHKERRQ(ierr);
        ierr = getEdgesFieldData(dataH1,meshPositionsFieldName); CHKERRQ(ierr);
        ierr = getTrisFieldData(dataH1,meshPositionsFieldName); CHKERRQ(ierr);
        ierr = getNodesFieldDofs(dataH1,meshPositionsFieldName); CHKERRQ(ierr);
        ierr = getEdgesFieldDofs(dataH1,meshPositionsFieldName); CHKERRQ(ierr);
        ierr = getTrisFieldDofs(dataH1,meshPositionsFieldName); CHKERRQ(ierr);
        try {
          ierr = opHONormals.opRhs(dataH1); CHKERRQ(ierr);
          ierr = opHONormals.calculateNormals(); CHKERRQ(ierr);
        } catch (exception& ex) {
          ostringstream ss;
          ss << "thorw in method: " << ex.what() << " at line " << __LINE__ << " in file " << __FILE__;
          SETERRQ(PETSC_COMM_SELF,1,ss.str().c_str());
        }
      }
    } catch (exception& ex) {
      ostringstream ss;
      ss << "thorw in method: " << ex.what() << " at line " << __LINE__ << " in file " << __FILE__;
      SETERRQ(PETSC_COMM_SELF,MOFEM_STD_EXCEPTION_THROW,ss.str().c_str());
    }

    if(dataH1.spacesOnEntities[MBTRI].test(HDIV)) {
      SETERRQ(PETSC_COMM_SELF,MOFEM_NOT_IMPLEMENTED,"not implemented yet");
    }

    /*for(int ss = 0;ss<2;ss++) {

      boost::ptr_vector<ForcesAndSurcesCore::UserDataOperator>::iterator oit,hi_oit;
      if(!ss) {

        oit = rowOpPtrVector.begin();
        hi_oit = rowOpPtrVector.end();

      } else {

        oit = colOpPtrVector.begin();
        hi_oit = colOpPtrVector.end();

      }

      for(;oit!=hi_oit;oit++) {

        oit->setPtrFE(this);
        BitFieldId data_id;
        FieldSpace space;
        string data_field_name;

        if(!ss) {
          data_id = mField.get_field_structure(oit->rowFieldName)->get_id();
          if((oit->getMoFEMFEPtr()->get_BitFieldId_data()&data_id).none()) {
            SETERRQ1(PETSC_COMM_SELF,MOFEM_DATA_INCONSISTENCT,"no data field < %s > on finite elemeny",oit->rowFieldName.c_str());
          }
          space = mField.get_field_structure(oit->rowFieldName)->get_space();
          data_field_name = oit->rowFieldName;
        } else {
          data_id = mField.get_field_structure(oit->colFieldName)->get_id();
          if((oit->getMoFEMFEPtr()->get_BitFieldId_data()&data_id).none()) {
            SETERRQ1(PETSC_COMM_SELF,MOFEM_DATA_INCONSISTENCT,"no data field < %s > on finite elemeny",oit->colFieldName.c_str());
          }
          space = mField.get_field_structure(oit->colFieldName)->get_space();
          data_field_name = oit->colFieldName;
        }

        DataForcesAndSurcesCore *op_data = NULL;
        switch(space) {
          case H1:
          op_data = &dataH1;
          break;
          case HCURL:
          SETERRQ(PETSC_COMM_SELF,MOFEM_NOT_IMPLEMENTED,"not implemented yet");
          break;
          case HDIV:
          op_data = &dataHdiv;
          break;
          case L2:
          SETERRQ(PETSC_COMM_SELF,MOFEM_NOT_IMPLEMENTED,"not implemented yet");
          break;
          case NOFIELD:
          op_data = &dataNoField;
          break;
          default:
          SETERRQ(PETSC_COMM_SELF,MOFEM_DATA_INCONSISTENCT,"data inconsistency");
          break;
        }

        try {

          switch(space) {
            case H1:
            if(!ss) {
              ierr = getRowNodesIndices(*op_data,data_field_name); CHKERRQ(ierr);
            } else {
              ierr = getColNodesIndices(*op_data,data_field_name); CHKERRQ(ierr);
            }
            ierr = getNodesFieldData(*op_data,data_field_name); CHKERRQ(ierr);
            ierr = getNodesFieldDofs(*op_data,data_field_name); CHKERRQ(ierr);
            case HCURL:
            if(!ss) {
              ierr = getEdgesRowIndices(*op_data,data_field_name); CHKERRQ(ierr);
            } else {
              ierr = getEdgesColIndices(*op_data,data_field_name); CHKERRQ(ierr);
            }
            ierr = getEdgesOrder(*op_data,data_field_name); CHKERRQ(ierr);
            ierr = getEdgesFieldData(*op_data,data_field_name); CHKERRQ(ierr);
            ierr = getEdgesFieldDofs(*op_data,data_field_name); CHKERRQ(ierr);
            case HDIV:
            case L2:
            if(!ss) {
              ierr = getTrisRowIndices(*op_data,data_field_name); CHKERRQ(ierr);
            } else {
              ierr = getTrisColIndices(*op_data,data_field_name); CHKERRQ(ierr);
            }
            ierr = getTrisOrder(*op_data,data_field_name); CHKERRQ(ierr);
            ierr = getTrisFieldData(*op_data,data_field_name); CHKERRQ(ierr);
            ierr = getTrisFieldDofs(*op_data,data_field_name); CHKERRQ(ierr);
            break;
            case NOFIELD:
            if(!getNinTheLoop()) {
              // NOFIELD data are the same for each element, can be retreived only once
              if(!ss) {
                ierr = getNoFieldRowIndices(*op_data,data_field_name); CHKERRQ(ierr);
              } else {
                ierr = getNoFieldColIndices(*op_data,data_field_name); CHKERRQ(ierr);
              }
              ierr = getNoFieldFieldData(*op_data,data_field_name); CHKERRQ(ierr);
              ierr = getNoFieldFieldDofs(*op_data,data_field_name); CHKERRQ(ierr);
            }
            break;
            default:
            break;
          }

        } catch (exception& ex) {
          ostringstream ss;
          ss << "thorw in method: " << ex.what() << " at line " << __LINE__ << " in file " << __FILE__;
          SETERRQ(PETSC_COMM_SELF,MOFEM_STD_EXCEPTION_THROW,ss.str().c_str());
        }

        try {
          ierr = oit->opRhs(*op_data); CHKERRQ(ierr);
        } catch (exception& ex) {
          ostringstream ss;
          ss << "thorw in method: " << ex.what() << " at line " << __LINE__ << " in file " << __FILE__;
          SETERRQ(PETSC_COMM_SELF,MOFEM_STD_EXCEPTION_THROW,ss.str().c_str());
        }

      }

    }

    for(
      boost::ptr_vector<ForcesAndSurcesCore::UserDataOperator>::iterator oit = rowColOpPtrVector.begin();
      oit != rowColOpPtrVector.end(); oit++
    ) {

      oit->setPtrFE(this);
      BitFieldId row_id = mField.get_field_structure(oit->rowFieldName)->get_id();
      BitFieldId col_id = mField.get_field_structure(oit->colFieldName)->get_id();

      if((oit->getMoFEMFEPtr()->get_BitFieldId_data()&row_id).none()) {
        SETERRQ1(PETSC_COMM_SELF,MOFEM_NOT_FOUND,"no row field < %s > on finite elemeny",oit->rowFieldName.c_str());
      }
      if((oit->getMoFEMFEPtr()->get_BitFieldId_data()&col_id).none()) {
        SETERRQ1(PETSC_COMM_SELF,MOFEM_NOT_FOUND,"no data field < %s > on finite elemeny",oit->colFieldName.c_str());
      }

      FieldSpace row_space = mField.get_field_structure(oit->rowFieldName)->get_space();

      DataForcesAndSurcesCore *row_op_data = NULL;
      switch(row_space) {
        case H1:
        row_op_data = &dataH1;
        break;
        case HCURL:
        SETERRQ(PETSC_COMM_SELF,MOFEM_NOT_IMPLEMENTED,"not implemented yet");
        break;
        case HDIV:
        row_op_data = &dataHdiv;
        break;
        case L2:
        SETERRQ(PETSC_COMM_SELF,MOFEM_NOT_IMPLEMENTED,"not implemented yet");
        break;
        case NOFIELD:
        row_op_data = &dataNoField;
        break;
        default:
        SETERRQ(PETSC_COMM_SELF,MOFEM_DATA_INCONSISTENCT,"data inconsistency");
        break;
      }

      switch(row_space) {
        case H1:
        ierr = getRowNodesIndices(*row_op_data,oit->rowFieldName); CHKERRQ(ierr);
        ierr = getNodesFieldData(*row_op_data,oit->rowFieldName); CHKERRQ(ierr);
        ierr = getNodesFieldDofs(*row_op_data,oit->rowFieldName); CHKERRQ(ierr);
        case HCURL:
        ierr = getEdgesRowIndices(*row_op_data,oit->rowFieldName); CHKERRQ(ierr);
        ierr = getEdgesOrder(*row_op_data,oit->rowFieldName); CHKERRQ(ierr);
        ierr = getEdgesFieldData(*row_op_data,oit->rowFieldName); CHKERRQ(ierr);
        ierr = getEdgesFieldDofs(*row_op_data,oit->rowFieldName); CHKERRQ(ierr);
        case HDIV:
        case L2:
        ierr = getTrisRowIndices(*row_op_data,oit->rowFieldName); CHKERRQ(ierr);
        ierr = getTrisOrder(*row_op_data,oit->rowFieldName); CHKERRQ(ierr);
        ierr = getTrisFieldData(*row_op_data,oit->rowFieldName); CHKERRQ(ierr);
        ierr = getTrisFieldDofs(*row_op_data,oit->rowFieldName); CHKERRQ(ierr);
        break;
        case NOFIELD:
        if(!getNinTheLoop()) {
          // NOFIELD data are the same for each element, can be retreived only once
          ierr = getNoFieldRowIndices(*row_op_data,oit->rowFieldName); CHKERRQ(ierr);
          ierr = getNoFieldFieldData(*row_op_data,oit->rowFieldName); CHKERRQ(ierr);
          ierr = getNoFieldFieldDofs(*row_op_data,oit->rowFieldName); CHKERRQ(ierr);
        }
        default:
        break;
      }

      FieldSpace col_space = mField.get_field_structure(oit->colFieldName)->get_space();
      DataForcesAndSurcesCore *col_op_data = NULL;
      switch(col_space) {
        case H1:
        col_op_data = &derivedDataH1;
        break;
        case HCURL:
        SETERRQ(PETSC_COMM_SELF,MOFEM_NOT_IMPLEMENTED,"not implemented yet");
        break;
        case HDIV:
        col_op_data = &derivedDataHdiv;
        break;
        case L2:
        SETERRQ(PETSC_COMM_SELF,MOFEM_NOT_IMPLEMENTED,"not implemented yet");
        break;
        default:
        SETERRQ(PETSC_COMM_SELF,MOFEM_DATA_INCONSISTENCT,"data inconsistency");
        break;
        case NOFIELD:
        col_op_data = &dataNoFieldCol;
        break;
      }

      switch(col_space) {
        case H1:
        ierr = getColNodesIndices(*col_op_data,oit->colFieldName); CHKERRQ(ierr);
        ierr = getNodesFieldData(*col_op_data,oit->colFieldName); CHKERRQ(ierr);
        ierr = getNodesFieldDofs(*col_op_data,oit->colFieldName); CHKERRQ(ierr);
        case HCURL:
        ierr = getEdgesColIndices(*col_op_data,oit->colFieldName); CHKERRQ(ierr);
        ierr = getEdgesOrder(*col_op_data,oit->colFieldName); CHKERRQ(ierr);
        ierr = getEdgesFieldData(*col_op_data,oit->colFieldName); CHKERRQ(ierr);
        ierr = getEdgesFieldDofs(*col_op_data,oit->colFieldName); CHKERRQ(ierr);
        case HDIV:
        case L2:
        ierr = getTrisColIndices(*col_op_data,oit->colFieldName); CHKERRQ(ierr);
        ierr = getTrisOrder(*col_op_data,oit->colFieldName); CHKERRQ(ierr);
        ierr = getTrisFieldData(*col_op_data,oit->colFieldName); CHKERRQ(ierr);
        ierr = getTrisFieldDofs(*col_op_data,oit->colFieldName); CHKERRQ(ierr);
        break;
        case NOFIELD:
        if(!getNinTheLoop()) {
          // NOFIELD data are the same for each element, can be retreived only once
          ierr = getNoFieldColIndices(*col_op_data,oit->colFieldName); CHKERRQ(ierr);
          ierr = getNoFieldFieldData(*col_op_data,oit->colFieldName); CHKERRQ(ierr);
          ierr = getNoFieldFieldDofs(*col_op_data,oit->colFieldName); CHKERRQ(ierr);
        }
        default:
        break;
      }

      try {
        ierr = oit->opLhs(*row_op_data,*col_op_data,oit->sYmm); CHKERRQ(ierr);
      } catch (exception& ex) {
        ostringstream ss;
        ss << "thorw in method: " << ex.what() << " at line " << __LINE__ << " in file " << __FILE__;
        SETERRQ(PETSC_COMM_SELF,MOFEM_STD_EXCEPTION_THROW,ss.str().c_str());
      }

    }*/

  } catch (exception& ex) {
    ostringstream ss;
    ss << "thorw in method: " << ex.what() << " at line " << __LINE__ << " in file " << __FILE__;
    SETERRQ(PETSC_COMM_SELF,MOFEM_STD_EXCEPTION_THROW,ss.str().c_str());
  }

  PetscFunctionReturn(0);
}

}
