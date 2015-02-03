/* Copyright (C) 2013, Lukasz Kaczmarczyk (likask AT wp.pl)
* --------------------------------------------------------------
*
* DESCRIPTION: FIXME
*
* This is not exactly procedure for linear elatic dynamics, since jacobian is
* evaluated at every time step and snes procedure is involved. However it is
* implemented like that, to test methodology for general nonlinear problem.
*
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

namespace MoFEM {

template<class T> 
void cOnstructor(DataForcesAndSurcesCore *data,EntityType type,T) {
    
    switch (type) {
      case MBTET:
	data->dataOnEntities[MBVERTEX].push_back(new T());
	for(int ee = 0;ee<6;ee++) {
	  data->dataOnEntities[MBEDGE].push_back(new T());
	}
	for(int ff = 0;ff<4;ff++) {
	  data->dataOnEntities[MBTRI].push_back(new T());
	}
	data->dataOnEntities[MBTET].push_back(new T());
	break;
      case MBTRI:
	data->dataOnEntities[MBVERTEX].push_back(new T());
	for(int ee = 0;ee<3;ee++) {
	  data->dataOnEntities[MBEDGE].push_back(new T());
	}
	data->dataOnEntities[MBTRI].push_back(new T());
	break;
      case MBEDGE:
	data->dataOnEntities[MBVERTEX].push_back(new T());
	data->dataOnEntities[MBEDGE].push_back(new T());
	break;
      case MBVERTEX:
	data->dataOnEntities[MBVERTEX].push_back(new T());
	break;
      case MBPRISM:
	data->dataOnEntities[MBVERTEX].push_back(new T());
	for(int ee = 0;ee<9;ee++) {
	  data->dataOnEntities[MBEDGE].push_back(new T());
	}
	for(int ff = 0;ff<5;ff++) {
	  data->dataOnEntities[MBTRI].push_back(new T());
	}
	break;
      default:
	throw("not implemenyed");
    }

}

DataForcesAndSurcesCore::DataForcesAndSurcesCore(EntityType type) {
  cOnstructor(this,type,EntData());
}


DerivedDataForcesAndSurcesCore::DerivedDataForcesAndSurcesCore(DataForcesAndSurcesCore &data): DataForcesAndSurcesCore() {

    boost::ptr_vector<EntData>::iterator iit;

    boost::ptr_vector<EntData>::iterator it;
    for(it = data.dataOnEntities[MBVERTEX].begin();it!=data.dataOnEntities[MBVERTEX].end();it++) {
      dataOnEntities[MBVERTEX].push_back(new DerivedEntData(*it));
    }
    for(it = data.dataOnEntities[MBEDGE].begin();it!=data.dataOnEntities[MBEDGE].end();it++) {
      dataOnEntities[MBEDGE].push_back(new DerivedEntData(*it));
    }
    for(it = data.dataOnEntities[MBTRI].begin();it!=data.dataOnEntities[MBTRI].end();it++) {
      dataOnEntities[MBTRI].push_back(new DerivedEntData(*it));
    }
    for(it = data.dataOnEntities[MBTET].begin();it!=data.dataOnEntities[MBTET].end();it++) {
      dataOnEntities[MBTET].push_back(new DerivedEntData(*it));
    }
  }

ostream& operator<<(ostream& os,const DataForcesAndSurcesCore::EntData &e) {
  os << 
    "sEnse: " << e.getSense() << endl << 
    "oRder: " << e.getOrder() << endl <<
    "iNdices: " << e.getIndices() << endl;
  os.precision(2);
  os << 
    "fieldData: " << std::fixed << e.getFieldData() << endl;
  os <<
    "N: " << std::fixed << e.getN() << endl <<
    "diffN: " << std::fixed << e.getDiffN();
  return os;
}

ostream& operator<<(ostream& os,const DataForcesAndSurcesCore &e) {
  for(unsigned int nn = 0;nn < e.dataOnEntities[MBVERTEX].size(); nn++) {
    os << "dataOnEntities[MBVERTEX][" << nn << "]" << endl << e.dataOnEntities[MBVERTEX][nn] << endl;
  }
  for(unsigned int ee = 0;ee < e.dataOnEntities[MBEDGE].size(); ee++) {
    os << "dataOnEntities[MBEDGE][" << ee << "]" << endl << e.dataOnEntities[MBEDGE][ee] << endl;
  }
  for(unsigned int ff = 0;ff < e.dataOnEntities[MBTRI].size(); ff++) {
    os << "dataOnEntities[MBTRI][" << ff << "] " << endl << e.dataOnEntities[MBTRI][ff] << endl;
  }
  for(unsigned int vv = 0;vv < e.dataOnEntities[MBTET].size(); vv++) {
    os << "dataOnEntities[MBTET][" << vv << "]" << endl << e.dataOnEntities[MBTET][vv] << endl;
  }
  return os;
}

PetscErrorCode ForcesAndSurcesCore::getSense(EntityType type,boost::ptr_vector<DataForcesAndSurcesCore::EntData> &data) {
    PetscFunctionBegin;
    SideNumber_multiIndex& side_table = const_cast<SideNumber_multiIndex&>(fePtr->get_side_number_table());
    if(data.size() < side_table.get<2>().count(type)) { // prims has 9 edges, some of edges for "flat" prism are not active
      SETERRQ2(PETSC_COMM_SELF,MOFEM_DATA_INSONSISTENCY,"data inconsistency %u < %u",data.size(),side_table.get<2>().count(type));
    }
    SideNumber_multiIndex::nth_index<2>::type::iterator siit = side_table.get<2>().lower_bound(type);
    SideNumber_multiIndex::nth_index<2>::type::iterator hi_siit = side_table.get<2>().upper_bound(type);
    for(;siit!=hi_siit;siit++) {
      data[siit->side_number].getSense() = siit->sense;
    }
    PetscFunctionReturn(0);
  }

PetscErrorCode ForcesAndSurcesCore::getEdgesSense(DataForcesAndSurcesCore &data) {
    PetscFunctionBegin;
    PetscErrorCode ierr;
    ierr = getSense(MBEDGE,data.dataOnEntities[MBEDGE]); CHKERRQ(ierr);
    PetscFunctionReturn(0);
  }

PetscErrorCode ForcesAndSurcesCore::getTrisSense(DataForcesAndSurcesCore &data) {
    PetscFunctionBegin;
    PetscErrorCode ierr;
    ierr = getSense(MBTRI,data.dataOnEntities[MBTRI]); CHKERRQ(ierr);
    PetscFunctionReturn(0);
  }

PetscErrorCode ForcesAndSurcesCore::getOrder(const EntityType type,const FieldSpace space,boost::ptr_vector<DataForcesAndSurcesCore::EntData> &data) {
    PetscFunctionBegin;
    SideNumber_multiIndex& side_table = const_cast<SideNumber_multiIndex&>(fePtr->get_side_number_table());
    if(data.size() < side_table.get<2>().count(type)) { // prims has 9 edges, some of edges for "flat" prism are not active
      SETERRQ2(PETSC_COMM_SELF,MOFEM_DATA_INSONSISTENCY,
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
	SETERRQ(PETSC_COMM_SELF,MOFEM_DATA_INSONSISTENCY,"data inconsistency");
      }
      data[side_number].getOrder() = data[side_number].getOrder() > ent_order ? data[side_number].getOrder() : ent_order;
    }
    PetscFunctionReturn(0);
  }


PetscErrorCode ForcesAndSurcesCore::getEdgesOrder(DataForcesAndSurcesCore &data,const FieldSpace space) {
    PetscFunctionBegin;
    PetscErrorCode ierr;
    ierr = getOrder(MBEDGE,space,data.dataOnEntities[MBEDGE]); CHKERRQ(ierr);
    PetscFunctionReturn(0);
  }

PetscErrorCode ForcesAndSurcesCore::getTrisOrder(DataForcesAndSurcesCore &data,const FieldSpace space) {
    PetscFunctionBegin;
    PetscErrorCode ierr;
    ierr = getOrder(MBTRI,space,data.dataOnEntities[MBTRI]); CHKERRQ(ierr);
    PetscFunctionReturn(0);
  }

PetscErrorCode ForcesAndSurcesCore::getTetsOrder(DataForcesAndSurcesCore &data,const FieldSpace space) {
    PetscFunctionBegin;
    PetscErrorCode ierr;
    ierr = getOrder(MBTET,space,data.dataOnEntities[MBTET]); CHKERRQ(ierr);
    PetscFunctionReturn(0);
  }

PetscErrorCode ForcesAndSurcesCore::getOrder(const string &field_name,const EntityType type,boost::ptr_vector<DataForcesAndSurcesCore::EntData> &data) {
    PetscFunctionBegin;
    SideNumber_multiIndex& side_table = const_cast<SideNumber_multiIndex&>(fePtr->get_side_number_table());
    if(data.size() < side_table.get<2>().count(type)) { // prims has 9 edges, some of edges for "flat" prism are not active
      SETERRQ2(PETSC_COMM_SELF,MOFEM_DATA_INSONSISTENCY,
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
	SETERRQ(PETSC_COMM_SELF,MOFEM_DATA_INSONSISTENCY,"data inconsistency");
      }
      data[side_number].getOrder() = data[side_number].getOrder() > ent_order ? data[side_number].getOrder() : ent_order;
    }
    PetscFunctionReturn(0);
  }

PetscErrorCode ForcesAndSurcesCore::getEdgesOrder(DataForcesAndSurcesCore &data,const string &field_name) {
    PetscFunctionBegin;
    PetscErrorCode ierr;
    ierr = getOrder(field_name,MBEDGE,data.dataOnEntities[MBEDGE]); CHKERRQ(ierr);
    PetscFunctionReturn(0);
  }

PetscErrorCode ForcesAndSurcesCore::getTrisOrder(DataForcesAndSurcesCore &data,const string &field_name) {
    PetscFunctionBegin;
    PetscErrorCode ierr;
    ierr = getOrder(field_name,MBTRI,data.dataOnEntities[MBTRI]); CHKERRQ(ierr);
    PetscFunctionReturn(0);
  }

PetscErrorCode ForcesAndSurcesCore::getTetsOrder(DataForcesAndSurcesCore &data,const string &field_name) {
    PetscFunctionBegin;
    PetscErrorCode ierr;
    ierr = getOrder(field_name,MBTET,data.dataOnEntities[MBTET]); CHKERRQ(ierr);
    PetscFunctionReturn(0);
  }

PetscErrorCode ForcesAndSurcesCore::getNodesIndices(
    const string &field_name,FENumeredDofMoFEMEntity_multiIndex &dofs,ublas::vector<int> &nodes_indices) {
    PetscFunctionBegin;
    FENumeredDofMoFEMEntity_multiIndex::index<Composite_Name_And_Type_mi_tag>::type::iterator dit,hi_dit;
    dit = dofs.get<Composite_Name_And_Type_mi_tag>().lower_bound(boost::make_tuple(field_name,MBVERTEX));
    hi_dit = dofs.get<Composite_Name_And_Type_mi_tag>().upper_bound(boost::make_tuple(field_name,MBVERTEX));
    nodes_indices.resize(distance(dit,hi_dit));
    for(;dit!=hi_dit;dit++) {
      int idx = dit->get_petsc_gloabl_dof_idx();
      int side_number = dit->side_number_ptr->side_number;
      nodes_indices[side_number*dit->get_max_rank()+dit->get_dof_rank()] = idx;
    }
    PetscFunctionReturn(0);
  }

PetscErrorCode ForcesAndSurcesCore::getRowNodesIndices(DataForcesAndSurcesCore &data,const string &field_name) {
    PetscFunctionBegin;
    PetscErrorCode ierr;
    ierr = getNodesIndices(
      field_name,const_cast<FENumeredDofMoFEMEntity_multiIndex&>(fePtr->get_rows_dofs()),data.dataOnEntities[MBVERTEX][0].getIndices()); CHKERRQ(ierr);
    PetscFunctionReturn(0);
  }

PetscErrorCode ForcesAndSurcesCore::getColNodesIndices(DataForcesAndSurcesCore &data,const string &field_name) {
    PetscFunctionBegin;
    PetscErrorCode ierr;
    ierr = getNodesIndices(field_name,const_cast<FENumeredDofMoFEMEntity_multiIndex&>(fePtr->get_cols_dofs()),data.dataOnEntities[MBVERTEX][0].getIndices()); CHKERRQ(ierr);
    PetscFunctionReturn(0);
  }

PetscErrorCode ForcesAndSurcesCore::getTypeIndices(
  const string &field_name,FENumeredDofMoFEMEntity_multiIndex &dofs,EntityType type,int side_number,ublas::vector<int> &indices) {
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
  boost::ptr_vector<DataForcesAndSurcesCore::EntData> &data) {
    PetscFunctionBegin;
    PetscErrorCode ierr;
    SideNumber_multiIndex& side_table = const_cast<SideNumber_multiIndex&>(fePtr->get_side_number_table());
    SideNumber_multiIndex::nth_index<2>::type::iterator siit = side_table.get<2>().lower_bound(type);
    SideNumber_multiIndex::nth_index<2>::type::iterator hi_siit = side_table.get<2>().upper_bound(type);
    for(;siit!=hi_siit;siit++) {
      ierr = getTypeIndices(field_name,dofs,type,siit->side_number,data[siit->side_number].getIndices()); CHKERRQ(ierr);
    }
    PetscFunctionReturn(0);
  }

PetscErrorCode ForcesAndSurcesCore::getEdgesRowIndices(DataForcesAndSurcesCore &data,const string &field_name) {
    PetscFunctionBegin;
    PetscErrorCode ierr;
    ierr = getTypeIndices(
      field_name,const_cast<FENumeredDofMoFEMEntity_multiIndex&>(fePtr->get_rows_dofs()),MBEDGE,data.dataOnEntities[MBEDGE]); CHKERRQ(ierr);
    PetscFunctionReturn(0);
  }

PetscErrorCode ForcesAndSurcesCore::getEdgesColIndices(DataForcesAndSurcesCore &data,const string &field_name) {
    PetscFunctionBegin;
    PetscErrorCode ierr;
    ierr = getTypeIndices(
      field_name,const_cast<FENumeredDofMoFEMEntity_multiIndex&>(fePtr->get_cols_dofs()),MBEDGE,data.dataOnEntities[MBEDGE]); CHKERRQ(ierr);
    PetscFunctionReturn(0);
  }

PetscErrorCode ForcesAndSurcesCore::getTrisRowIndices(
    DataForcesAndSurcesCore &data,const string &field_name) {
    PetscFunctionBegin;
    PetscErrorCode ierr;
    ierr = getTypeIndices(field_name,const_cast<FENumeredDofMoFEMEntity_multiIndex&>(fePtr->get_rows_dofs()),MBTRI,data.dataOnEntities[MBTRI]); CHKERRQ(ierr);
    PetscFunctionReturn(0);
  }

PetscErrorCode ForcesAndSurcesCore::getTrisColIndices(DataForcesAndSurcesCore &data,const string &field_name) {
    PetscFunctionBegin;
    PetscErrorCode ierr;
    ierr = getTypeIndices(field_name,const_cast<FENumeredDofMoFEMEntity_multiIndex&>(fePtr->get_cols_dofs()),MBTRI,data.dataOnEntities[MBTRI]); CHKERRQ(ierr);
    PetscFunctionReturn(0);
  }

PetscErrorCode ForcesAndSurcesCore::getTetsRowIndices(DataForcesAndSurcesCore &data,const string &field_name) {
    PetscFunctionBegin;
    PetscErrorCode ierr;
    if(data.dataOnEntities[MBTET].size() == 0) {
      SETERRQ(PETSC_COMM_SELF,MOFEM_DATA_INSONSISTENCY,"data inconsistency");
    }
    ierr = getTypeIndices(field_name,const_cast<FENumeredDofMoFEMEntity_multiIndex&>(fePtr->get_rows_dofs()),MBTET,0,data.dataOnEntities[MBTET][0].getIndices()); CHKERRQ(ierr);
    PetscFunctionReturn(0);
  }

PetscErrorCode ForcesAndSurcesCore::getTetsColIndices(DataForcesAndSurcesCore &data,const string &field_name) {
    PetscFunctionBegin;
    PetscErrorCode ierr;
    if(data.dataOnEntities[MBTET].size() == 0) {
      SETERRQ(PETSC_COMM_SELF,MOFEM_DATA_INSONSISTENCY,"data inconsistency");
    }
    ierr = getTypeIndices(field_name,const_cast<FENumeredDofMoFEMEntity_multiIndex&>(fePtr->get_cols_dofs()),MBTET,0,data.dataOnEntities[MBTET][0].getIndices()); CHKERRQ(ierr);
    PetscFunctionReturn(0);
  }

PetscErrorCode ForcesAndSurcesCore::getNodesFieldData(
    const string &field_name,FEDofMoFEMEntity_multiIndex &dofs,ublas::vector<FieldData> &nodes_field_data) {
    PetscFunctionBegin;
    FEDofMoFEMEntity_multiIndex::index<Composite_Name_And_Type_mi_tag>::type::iterator dit,hi_dit;
    dit = dofs.get<Composite_Name_And_Type_mi_tag>().lower_bound(boost::make_tuple(field_name,MBVERTEX));
    hi_dit = dofs.get<Composite_Name_And_Type_mi_tag>().upper_bound(boost::make_tuple(field_name,MBVERTEX));
    nodes_field_data.resize(distance(dit,hi_dit));
    for(;dit!=hi_dit;dit++) {
      FieldData val = dit->get_FieldData();
      int side_number = dit->side_number_ptr->side_number;
      nodes_field_data[side_number*dit->get_max_rank()+dit->get_dof_rank()] = val;
    }
    PetscFunctionReturn(0);
  }

PetscErrorCode ForcesAndSurcesCore::getNodesFieldData(DataForcesAndSurcesCore &data,const string &field_name) {
    PetscFunctionBegin;
    PetscErrorCode ierr;
    ierr = getNodesFieldData(
      field_name,const_cast<FEDofMoFEMEntity_multiIndex&>(fePtr->get_data_dofs()),data.dataOnEntities[MBVERTEX][0].getFieldData()); CHKERRQ(ierr);
    PetscFunctionReturn(0);
  }

PetscErrorCode ForcesAndSurcesCore::getTypeFieldData(
  const string &field_name,FEDofMoFEMEntity_multiIndex &dofs,EntityType type,int side_number,ublas::vector<FieldData> &ent_field_data) {
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
  boost::ptr_vector<DataForcesAndSurcesCore::EntData> &data) {
    PetscFunctionBegin;
    PetscErrorCode ierr;
    SideNumber_multiIndex& side_table = const_cast<SideNumber_multiIndex&>(fePtr->get_side_number_table());
    if(data.size() < side_table.get<2>().count(type)) { 
      SETERRQ(PETSC_COMM_SELF,MOFEM_DATA_INSONSISTENCY,"data inconsistency");
    }
    SideNumber_multiIndex::nth_index<2>::type::iterator siit = side_table.get<2>().lower_bound(type);
    SideNumber_multiIndex::nth_index<2>::type::iterator hi_siit = side_table.get<2>().upper_bound(type);
    for(;siit!=hi_siit;siit++) {
      ierr = getTypeFieldData(field_name,dofs,type,siit->side_number,data[siit->side_number].getFieldData()); CHKERRQ(ierr);
    }
    PetscFunctionReturn(0);
  }

PetscErrorCode ForcesAndSurcesCore::getEdgesFieldData(DataForcesAndSurcesCore &data,const string &field_name) {
    PetscFunctionBegin;
    PetscErrorCode ierr;
    ierr = getTypeFieldData(
      field_name,const_cast<FEDofMoFEMEntity_multiIndex&>(fePtr->get_data_dofs()),MBEDGE,data.dataOnEntities[MBEDGE]); CHKERRQ(ierr);
    PetscFunctionReturn(0);
  }

PetscErrorCode ForcesAndSurcesCore::getTrisFieldData(
    DataForcesAndSurcesCore &data,const string &field_name) {
    PetscFunctionBegin;
    PetscErrorCode ierr;
    ierr = getTypeFieldData(field_name,const_cast<FEDofMoFEMEntity_multiIndex&>(fePtr->get_data_dofs()),MBTRI,data.dataOnEntities[MBTRI]); CHKERRQ(ierr);
    PetscFunctionReturn(0);
  }

PetscErrorCode ForcesAndSurcesCore::getTetsFieldData(DataForcesAndSurcesCore &data,const string &field_name) {
    PetscFunctionBegin;
    PetscErrorCode ierr;
    if(data.dataOnEntities[MBTET].size() == 0) {
      SETERRQ(PETSC_COMM_SELF,MOFEM_DATA_INSONSISTENCY,"data inconsistency");
    }
    ierr = getTypeFieldData(field_name,const_cast<FEDofMoFEMEntity_multiIndex&>(fePtr->get_data_dofs()),MBTET,0,data.dataOnEntities[MBTET][0].getFieldData()); CHKERRQ(ierr);
    PetscFunctionReturn(0);
  }

PetscErrorCode ForcesAndSurcesCore::getNodesFieldDofs(
    const string &field_name,FEDofMoFEMEntity_multiIndex &dofs,ublas::vector<const FEDofMoFEMEntity*> &nodes_field_dofs) {
    PetscFunctionBegin;
    FEDofMoFEMEntity_multiIndex::index<Composite_Name_And_Type_mi_tag>::type::iterator dit,hi_dit;
    dit = dofs.get<Composite_Name_And_Type_mi_tag>().lower_bound(boost::make_tuple(field_name,MBVERTEX));
    hi_dit = dofs.get<Composite_Name_And_Type_mi_tag>().upper_bound(boost::make_tuple(field_name,MBVERTEX));
    nodes_field_dofs.resize(distance(dit,hi_dit));
    for(;dit!=hi_dit;dit++) {
      int side_number = dit->side_number_ptr->side_number;
      nodes_field_dofs[side_number*dit->get_max_rank()+dit->get_dof_rank()] = &*dit;
    }
    PetscFunctionReturn(0);
  }

PetscErrorCode ForcesAndSurcesCore::getNodesFieldDofs(DataForcesAndSurcesCore &data,const string &field_name) {
    PetscFunctionBegin;
    PetscErrorCode ierr;
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
    PetscErrorCode ierr;
    SideNumber_multiIndex& side_table = const_cast<SideNumber_multiIndex&>(fePtr->get_side_number_table());
    if(data.size() < side_table.get<2>().count(type)) {
      SETERRQ(PETSC_COMM_SELF,MOFEM_DATA_INSONSISTENCY,"data inconsistency");
    }
    SideNumber_multiIndex::nth_index<2>::type::iterator siit = side_table.get<2>().lower_bound(type);
    SideNumber_multiIndex::nth_index<2>::type::iterator hi_siit = side_table.get<2>().upper_bound(type);
    for(;siit!=hi_siit;siit++) {
      ierr = getTypeFieldDofs(field_name,dofs,type,siit->side_number,data[siit->side_number].getFieldDofs()); CHKERRQ(ierr);
    }
    PetscFunctionReturn(0);
  }

PetscErrorCode ForcesAndSurcesCore::getEdgesFieldDofs(DataForcesAndSurcesCore &data,const string &field_name) {
    PetscFunctionBegin;
    PetscErrorCode ierr;
    ierr = getTypeFieldDofs(
      field_name,const_cast<FEDofMoFEMEntity_multiIndex&>(fePtr->get_data_dofs()),MBEDGE,data.dataOnEntities[MBEDGE]); CHKERRQ(ierr);
    PetscFunctionReturn(0);
  }

PetscErrorCode ForcesAndSurcesCore::getTrisFieldDofs(DataForcesAndSurcesCore &data,const string &field_name) {
    PetscFunctionBegin;
    PetscErrorCode ierr;
    ierr = getTypeFieldDofs(field_name,const_cast<FEDofMoFEMEntity_multiIndex&>(fePtr->get_data_dofs()),MBTRI,data.dataOnEntities[MBTRI]); CHKERRQ(ierr);
    PetscFunctionReturn(0);
  }

PetscErrorCode ForcesAndSurcesCore::getTetsFieldDofs(DataForcesAndSurcesCore &data,const string &field_name) {
    PetscFunctionBegin;
    PetscErrorCode ierr;
    if(data.dataOnEntities[MBTET].size() == 0) {
      SETERRQ(PETSC_COMM_SELF,MOFEM_DATA_INSONSISTENCY,"data inconsistency");
    }
    ierr = getTypeFieldDofs(field_name,const_cast<FEDofMoFEMEntity_multiIndex&>(fePtr->get_data_dofs()),MBTET,0,data.dataOnEntities[MBTET][0].getFieldDofs()); CHKERRQ(ierr);
    PetscFunctionReturn(0);
  }

PetscErrorCode ForcesAndSurcesCore::getFaceNodes(DataForcesAndSurcesCore &data) {
    PetscFunctionBegin;
    //PetscAttachDebugger();
    ErrorCode rval;
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
	if(num_nodes_tet != 4) SETERRQ(PETSC_COMM_SELF,MOFEM_DATA_INSONSISTENCY,"data inconsistency");
	int num_nodes_face;
	const EntityHandle *conn_face;
	rval = mField.get_moab().get_connectivity(side->ent,conn_face,num_nodes_face,true); CHKERR_PETSC(rval);
	if(num_nodes_face != 3) SETERRQ(PETSC_COMM_SELF,1,"data inconsistency");
	if(conn_face[0] != conn_tet[data.facesNodes(side->side_number,0)]) SETERRQ(PETSC_COMM_SELF,MOFEM_DATA_INSONSISTENCY,"data inconsistency");
	if(conn_face[1] != conn_tet[data.facesNodes(side->side_number,1)]) SETERRQ(PETSC_COMM_SELF,MOFEM_DATA_INSONSISTENCY,"data inconsistency");
	if(conn_face[2] != conn_tet[data.facesNodes(side->side_number,2)]) SETERRQ(PETSC_COMM_SELF,MOFEM_DATA_INSONSISTENCY,"data inconsistency");
      }
    }
    PetscFunctionReturn(0);
  }

PetscErrorCode ForcesAndSurcesCore::getSpacesOnEntities(DataForcesAndSurcesCore &data) {
  PetscFunctionBegin;

  for(_IT_GET_FEDATA_DOFS_FOR_LOOP_(this,dof)) {
    data.spacesOnEntities[dof->get_ent_type()].set(dof->get_space());
  }

  PetscFunctionReturn(0);
}

PetscErrorCode ForcesAndSurcesCore::shapeTETFunctions_H1(
    DataForcesAndSurcesCore &data,
    const double *G_X,const double *G_Y,const double *G_Z,const int G_DIM) {
    PetscFunctionBegin;

    PetscErrorCode ierr;

    data.dataOnEntities[MBVERTEX][0].getN().resize(G_DIM,4);
    ierr = ShapeMBTET(&*data.dataOnEntities[MBVERTEX][0].getN().data().begin(),G_X,G_Y,G_Z,G_DIM); CHKERRQ(ierr);
    data.dataOnEntities[MBVERTEX][0].getDiffN().resize(4,3);
    ierr = ShapeDiffMBTET(&*data.dataOnEntities[MBVERTEX][0].getDiffN().data().begin()); CHKERRQ(ierr);

    if((data.spacesOnEntities[MBEDGE]).test(H1)) {

    //edges
    if(data.dataOnEntities[MBEDGE].size()!=6) {
      SETERRQ(PETSC_COMM_SELF,MOFEM_DATA_INSONSISTENCY,"data inconsistency");
    }
    int _sense_[6],_order_[6];
    double *_H1edgeN_[6],*_diffH1edgeN_[6];
    for(int ee = 0;ee<6;ee++) {
      if(data.dataOnEntities[MBEDGE][ee].getSense() == 0) {
	SETERRQ(PETSC_COMM_SELF,MOFEM_DATA_INSONSISTENCY,"data inconsistency");
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
      SETERRQ(PETSC_COMM_SELF,MOFEM_DATA_INSONSISTENCY,"data inconsistency");
    }
    double *_H1faceN_[4],*_diffH1faceN_[4];
    for(int ff = 0;ff<4;ff++) {
      if(data.dataOnEntities[MBTRI][ff].getSense() == 0) {
	SETERRQ(PETSC_COMM_SELF,MOFEM_DATA_INSONSISTENCY,"data inconsistency");
      }
      int nb_dofs = NBFACE_H1(data.dataOnEntities[MBTRI][ff].getOrder());
      _order_[ff] = data.dataOnEntities[MBTRI][ff].getOrder();
      data.dataOnEntities[MBTRI][ff].getN().resize(G_DIM,nb_dofs);
      data.dataOnEntities[MBTRI][ff].getDiffN().resize(G_DIM,3*nb_dofs);
      _H1faceN_[ff] = &*data.dataOnEntities[MBTRI][ff].getN().data().begin();
      _diffH1faceN_[ff] = &*data.dataOnEntities[MBTRI][ff].getDiffN().data().begin();
    }
    if(data.facesNodes.size1() != 4) {
      SETERRQ(PETSC_COMM_SELF,MOFEM_DATA_INSONSISTENCY,"data inconsistency");
    }
    if(data.facesNodes.size2() != 3) {
      SETERRQ(PETSC_COMM_SELF,MOFEM_DATA_INSONSISTENCY,"data inconsistency");
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
      &*data.dataOnEntities[MBTET][0].getN().data().begin(),&*data.dataOnEntities[MBTET][0].getDiffN().data().begin(),G_DIM); CHKERRQ(ierr);

    }

    PetscFunctionReturn(0);
  }
  
PetscErrorCode ForcesAndSurcesCore::shapeTETFunctions_L2(
    DataForcesAndSurcesCore &data,
    const double *G_X,const double *G_Y,const double *G_Z,const int G_DIM) {
    PetscFunctionBegin;

    PetscErrorCode ierr;

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

    PetscErrorCode ierr;

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
	SETERRQ(PETSC_COMM_SELF,MOFEM_DATA_INSONSISTENCY,"data inconsistency");
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
      SETERRQ(PETSC_COMM_SELF,MOFEM_DATA_INSONSISTENCY,"data inconsistency");
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
  PetscErrorCode ierr;

  data.dataOnEntities[MBVERTEX][0].getN().resize(G_DIM,3);
  ierr = ShapeMBTRI(&*data.dataOnEntities[MBVERTEX][0].getN().data().begin(),G_X,G_Y,G_DIM); CHKERRQ(ierr);
  data.dataOnEntities[MBVERTEX][0].getDiffN().resize(3,2);
  ierr = ShapeDiffMBTRI(&*data.dataOnEntities[MBVERTEX][0].getDiffN().data().begin()); CHKERRQ(ierr);

  if((data.spacesOnEntities[MBEDGE]).test(H1)) {

  //edges
  if(data.dataOnEntities[MBEDGE].size()!=3) {
    SETERRQ(PETSC_COMM_SELF,MOFEM_DATA_INSONSISTENCY,"data inconsistency");
  }
  int _sense_[3],_order_[3];
  double *_H1edgeN_[3],*_diffH1edgeN_[3];
  for(int ee = 0;ee<3;ee++) {
    if(data.dataOnEntities[MBEDGE][ee].getSense() == 0) {
	SETERRQ(PETSC_COMM_SELF,MOFEM_DATA_INSONSISTENCY,"data inconsistency");
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
    SETERRQ(PETSC_COMM_SELF,MOFEM_DATA_INSONSISTENCY,"data inconsistency");
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

  PetscErrorCode ierr;

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
    SETERRQ(PETSC_COMM_SELF,MOFEM_DATA_INSONSISTENCY,"data inconsistency");
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
  PetscErrorCode ierr;

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
  PetscErrorCode ierr;

  data.dataOnEntities[MBVERTEX][0].getN().resize(G_DIM,6);
  ierr = ShapeMBTRI(&*data.dataOnEntities[MBVERTEX][0].getN().data().begin(),G_X,G_Y,G_DIM); CHKERRQ(ierr);
  data.dataOnEntities[MBVERTEX][0].getDiffN().resize(6,2);
  ierr = ShapeDiffMBTRI(&*data.dataOnEntities[MBVERTEX][0].getDiffN().data().begin()); CHKERRQ(ierr);
  //shape functions on other side are like on the first one
  for(int nn = 0;nn<3;nn++) {
    for(int gg = 0;gg<3;gg++) {
      data.dataOnEntities[MBVERTEX][0].getN()(gg,nn+3) = data.dataOnEntities[MBVERTEX][0].getN()(gg,nn);
    }
    for(int dd = 0;dd<2;dd++) {
      data.dataOnEntities[MBVERTEX][0].getDiffN()(nn+3,dd) = data.dataOnEntities[MBVERTEX][0].getDiffN()(nn,dd);
    }
  }

  //edges
  if(data.dataOnEntities[MBEDGE].size()!=9) {
    SETERRQ(PETSC_COMM_SELF,MOFEM_DATA_INSONSISTENCY,"data inconsistency");
  }
  
  int valid_edges[] = { 1,1,1, 0,0,0, 1,1,1 };
  int sense[9],order[9];
  double *H1edgeN[9],*diffH1edgeN[9];
  for(int ee = 0;ee<9;ee++) {
    if(!valid_edges[ee]) continue;
    if(data.dataOnEntities[MBEDGE][ee].getSense() == 0) {
	SETERRQ(PETSC_COMM_SELF,MOFEM_DATA_INSONSISTENCY,"data inconsistency");
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
    &*data.dataOnEntities[MBVERTEX][0].getN().data().begin(),&*data.dataOnEntities[MBVERTEX][0].getDiffN().data().begin(),
    &H1edgeN[0],&diffH1edgeN[0],G_DIM); CHKERRQ(ierr);
  //shape functions on face 4
  ierr = H1_EdgeShapeFunctions_MBTRI(&sense[6],&order[6],
    &*data.dataOnEntities[MBVERTEX][0].getN().data().begin(),&*data.dataOnEntities[MBVERTEX][0].getDiffN().data().begin(),
    &H1edgeN[6],&diffH1edgeN[6],G_DIM); CHKERRQ(ierr);

  //face
  if(data.dataOnEntities[MBTRI].size()!=5) {
    SETERRQ(PETSC_COMM_SELF,MOFEM_DATA_INSONSISTENCY,"data inconsistency");
  }
  for(int ff = 3;ff<5;ff++) {
    int nb_dofs = NBFACE_H1(data.dataOnEntities[MBTRI][ff].getOrder());
    data.dataOnEntities[MBTRI][ff].getN().resize(G_DIM,nb_dofs);
    data.dataOnEntities[MBTRI][ff].getDiffN().resize(G_DIM,2*nb_dofs);
    const int face_nodes[] = { 0,1,2 };
    ierr = H1_FaceShapeFunctions_MBTRI(
      face_nodes,data.dataOnEntities[MBTRI][ff].getOrder(),
      &*data.dataOnEntities[MBVERTEX][0].getN().data().begin(),&*data.dataOnEntities[MBVERTEX][0].getDiffN().data().begin(),
      &*data.dataOnEntities[MBTRI][ff].getN().data().begin(),&*data.dataOnEntities[MBTRI][ff].getDiffN().data().begin(),
      G_DIM); CHKERRQ(ierr);
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

}
