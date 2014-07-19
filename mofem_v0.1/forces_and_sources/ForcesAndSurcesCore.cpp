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


#include "ForcesAndSurcesCore.hpp"

#include<FEM.h>
#include<H1HdivHcurlL2.h>

extern "C" {
#include "gm_rule.h"
}

using namespace boost::numeric;

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
    if(data.size() != side_table.get<2>().count(type)) {
      SETERRQ2(PETSC_COMM_SELF,MOFEM_DATA_INSONSISTENCY,"data inconsistency %u != %u",data.size(),side_table.get<2>().count(type));
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
    if(data.size() != side_table.get<2>().count(type)) {
      SETERRQ2(PETSC_COMM_SELF,MOFEM_DATA_INSONSISTENCY,
	"data inconsistency %d != %d",
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
    if(data.size() != side_table.get<2>().count(type)) {
      SETERRQ2(PETSC_COMM_SELF,MOFEM_DATA_INSONSISTENCY,
	"data inconsistency %d != %d",
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
    if(data.size() != side_table.get<2>().count(type)) {
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

    //calulate shape function for tet, needed to construct shape functions for h_div 

    if(data.dataOnEntities[MBVERTEX][0].getN().size1() != (unsigned int)G_DIM) {
      data.dataOnEntities[MBVERTEX][0].getN().resize(G_DIM,4);
      ierr = ShapeMBTET(&*data.dataOnEntities[MBVERTEX][0].getN().data().begin(),G_X,G_Y,G_Z,G_DIM); CHKERRQ(ierr);
    }
    //that is cheep to calate, no harm done if recalulated
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


PetscErrorCode DataOperator::opLhs(
    DataForcesAndSurcesCore &row_data,DataForcesAndSurcesCore &col_data,bool symm) {
  PetscFunctionBegin;
  PetscErrorCode ierr;

  //nodes
  ierr = doWork(
    -1,-1,MBVERTEX,MBVERTEX,
    row_data.dataOnEntities[MBVERTEX][0],col_data.dataOnEntities[MBVERTEX][0]); CHKERRQ(ierr);
  for(unsigned int VV = 0;VV<col_data.dataOnEntities[MBTET].size();VV++) {
    ierr = doWork(
      -1,VV,MBVERTEX,MBTET,
      row_data.dataOnEntities[MBVERTEX][0],col_data.dataOnEntities[MBTET][VV]); CHKERRQ(ierr);
    if(!symm) {
      for(unsigned int EE = 0;EE<col_data.dataOnEntities[MBEDGE].size();EE++) {
	ierr = doWork(
	  -1,EE,MBVERTEX,MBEDGE,
	  row_data.dataOnEntities[MBVERTEX][0],col_data.dataOnEntities[MBEDGE][0]); CHKERRQ(ierr);
      }
      for(unsigned int FF = 0;FF<col_data.dataOnEntities[MBTRI].size();FF++) {
	ierr = doWork(
	  -1,FF,MBVERTEX,MBTRI,
	  row_data.dataOnEntities[MBVERTEX][0],col_data.dataOnEntities[MBTRI][FF]); CHKERRQ(ierr);
      }
    }
  }

  //edges
  for(unsigned int ee = 0;ee<row_data.dataOnEntities[MBEDGE].size();ee++) {
    if(row_data.dataOnEntities[MBEDGE][ee].getN().size1()==0) continue;
    ierr = doWork(
	ee,-1,MBEDGE,MBVERTEX,
	row_data.dataOnEntities[MBEDGE][ee],col_data.dataOnEntities[MBVERTEX][0]); CHKERRQ(ierr);
    for(unsigned int VV = 0;VV<col_data.dataOnEntities[MBTET].size();VV++) {
      ierr = doWork(
	ee,VV,MBEDGE,MBTET,
	row_data.dataOnEntities[MBEDGE][ee],col_data.dataOnEntities[MBTET][VV]); CHKERRQ(ierr);
    }
    unsigned int EE = 0;
    if(symm) EE = ee;
    for(;EE<col_data.dataOnEntities[MBEDGE].size();EE++) {
      ierr = doWork(
	  ee,EE,MBEDGE,MBEDGE,
	  row_data.dataOnEntities[MBEDGE][ee],col_data.dataOnEntities[MBEDGE][EE]); CHKERRQ(ierr);
    }
    for(unsigned int FF = 0;FF<col_data.dataOnEntities[MBTRI].size();FF++) {
      ierr = doWork(
	  ee,FF,MBEDGE,MBTRI,
	  row_data.dataOnEntities[MBEDGE][ee],col_data.dataOnEntities[MBTRI][FF]); CHKERRQ(ierr);
    }
  }

  //faces
  for(unsigned int ff = 0;ff<row_data.dataOnEntities[MBTRI].size();ff++) {
    ierr = doWork(
	ff,-1,MBTRI,MBVERTEX,
	row_data.dataOnEntities[MBTRI][ff],col_data.dataOnEntities[MBVERTEX][0]); CHKERRQ(ierr);
    for(unsigned int VV = 0;VV<col_data.dataOnEntities[MBTET].size();VV++) {
      ierr = doWork(
	ff,VV,MBTRI,MBTET,
	row_data.dataOnEntities[MBTRI][ff],col_data.dataOnEntities[MBTET][VV]); CHKERRQ(ierr);
    }
    unsigned int FF = 0;
    if(symm) FF = ff;
    for(;FF<col_data.dataOnEntities[MBTRI].size();FF++) {
      ierr = doWork(
	ff,FF,MBTRI,MBTRI,
	row_data.dataOnEntities[MBTRI][ff],col_data.dataOnEntities[MBTRI][FF]); CHKERRQ(ierr);
    }
    if(!symm) {
      unsigned int EE = 0;
      for(;EE<col_data.dataOnEntities[MBEDGE].size();EE++) {
	ierr = doWork(
	  ff,EE,MBTRI,MBEDGE,
	  row_data.dataOnEntities[MBTRI][ff],col_data.dataOnEntities[MBEDGE][EE]); CHKERRQ(ierr);
      }
    }
  }

  //volumes
  for(unsigned int vv = 0;vv<row_data.dataOnEntities[MBTET].size();vv++) {
    unsigned int VV = 0;
    if(symm) VV = vv;
    for(;VV<col_data.dataOnEntities[MBTET].size();VV++) {
      ierr = doWork(
	vv,VV,MBTET,MBTET,
	row_data.dataOnEntities[MBTET][vv],col_data.dataOnEntities[MBTET][VV]); CHKERRQ(ierr);
    }
    if(!symm) {
      //vertex
      ierr = doWork(
	VV,-1,MBTET,MBVERTEX,
	row_data.dataOnEntities[MBTET][vv],col_data.dataOnEntities[MBVERTEX][0]); CHKERRQ(ierr);
      //edges
      for(unsigned int EE = 0;EE<col_data.dataOnEntities[MBEDGE].size();EE++) {
	ierr = doWork(
	  VV,EE,MBTET,MBEDGE,
	  row_data.dataOnEntities[MBTET][vv],col_data.dataOnEntities[MBEDGE][EE]); CHKERRQ(ierr);
      }
      //faces
      for(unsigned int FF = 0;FF<col_data.dataOnEntities[MBTRI].size();FF++) {
	ierr = doWork(
	  VV,FF,MBTET,MBTRI,
	  row_data.dataOnEntities[MBTET][vv],col_data.dataOnEntities[MBTRI][FF]); CHKERRQ(ierr);
      }
    }
  }

  PetscFunctionReturn(0);
}

PetscErrorCode DataOperator::opRhs(DataForcesAndSurcesCore &data) {
  PetscFunctionBegin;
  PetscErrorCode ierr;

  for(unsigned int nn = 0;nn<data.dataOnEntities[MBVERTEX].size();nn++) {
    ierr = doWork(-1,MBVERTEX,data.dataOnEntities[MBVERTEX][nn]); CHKERRQ(ierr);
  }
  for(unsigned int ee = 0;ee<data.dataOnEntities[MBEDGE].size();ee++) {
    ierr = doWork(ee,MBEDGE,data.dataOnEntities[MBEDGE][ee]); CHKERRQ(ierr);
  }
  for(unsigned int ff = 0;ff<data.dataOnEntities[MBTRI].size();ff++) {
    ierr = doWork(ff,MBTRI,data.dataOnEntities[MBTRI][ff]); CHKERRQ(ierr);
  }
  for(unsigned int vv = 0;vv<data.dataOnEntities[MBTET].size();vv++) {
    ierr = doWork(-1,MBTET,data.dataOnEntities[MBTET][vv]); CHKERRQ(ierr);
  }

  PetscFunctionReturn(0);
}

PetscErrorCode OpSetInvJacH1::doWork(
    int side,
    EntityType type,
    DataForcesAndSurcesCore::EntData &data) {
  PetscFunctionBegin;
  PetscErrorCode ierr;

  try {

    diffNinvJac.resize(data.getDiffN().size1(),data.getDiffN().size2());
    unsigned int nb_gauss_pts = data.getN().size1();
    unsigned int nb_dofs = data.getN().size2();
    if(type!=MBVERTEX) {
      if(nb_dofs != data.getDiffN().size2()/3) {
        SETERRQ2(PETSC_COMM_SELF,MOFEM_DATA_INSONSISTENCY,
  	"data inconsistency nb_dofs != data.diffN.size2()/3 ( %u != %u/3 )",
  	nb_dofs,data.getDiffN().size2());
      }
    }
  
    switch (type) {
  
      case MBVERTEX: {
        ierr = ShapeDiffMBTETinvJ(
  	&*data.getDiffN().data().begin(),&*invJac.data().begin(),&*diffNinvJac.data().begin()); CHKERRQ(ierr);
      }
      break;
      case MBEDGE:
      case MBTRI:
      case MBTET: {
        for(unsigned int gg = 0;gg<nb_gauss_pts;gg++) {
	  for(unsigned int dd = 0;dd<nb_dofs;dd++) {
	    cblas_dgemv(CblasRowMajor,CblasTrans,3,3,1.,
	      &*invJac.data().begin(),3,&data.getDiffN()(gg,3*dd),1,0.,&diffNinvJac(gg,3*dd),1); 
	  }
        }
      }
      break;
      default:
        SETERRQ(PETSC_COMM_SELF,MOFEM_NOT_IMPLEMENTED,"not implemented");
  
    }
  
    data.getDiffN().data().swap(diffNinvJac.data());

  } catch (exception& ex) {
    ostringstream ss;
    ss << "thorw in method: " << ex.what() << " at line " << __LINE__ << " in file " << __FILE__;
    SETERRQ(PETSC_COMM_SELF,MOFEM_STD_EXCEPTION_THROW,ss.str().c_str());
  }

  PetscFunctionReturn(0);
}

PetscErrorCode OpSetInvJacHdiv::doWork(
    int side,
    EntityType type,
    DataForcesAndSurcesCore::EntData &data) {
  PetscFunctionBegin;

  if(type != MBTRI && type != MBTET) PetscFunctionReturn(0);

  try {

    diffHdiv_invJac.resize(data.getDiffHdivN().size1(),data.getDiffHdivN().size2());

    unsigned int nb_gauss_pts = data.getDiffHdivN().size1();
    unsigned int nb_dofs = data.getDiffHdivN().size2()/9;
    
    unsigned int gg = 0;
    for(;gg<nb_gauss_pts;gg++) {
      unsigned int dd = 0;
      for(;dd<nb_dofs;dd++) {
	const double *DiffHdivN = &((data.getDiffHdivN(gg))(dd,0));
	for(int kk = 0;kk<3;kk++) {
	  cblas_dgemv(CblasRowMajor,CblasTrans,3,3,1.,
	    &*invJac.data().begin(),3,&DiffHdivN[kk],3,
	    0.,&diffHdiv_invJac(gg,9*dd+kk),3); 
	}
      }
    }

    data.getDiffHdivN().data().swap(diffHdiv_invJac.data());

  } catch (exception& ex) {
    ostringstream ss;
    ss << "thorw in method: " << ex.what() << " at line " << __LINE__ << " in file " << __FILE__;
    SETERRQ(PETSC_COMM_SELF,MOFEM_STD_EXCEPTION_THROW,ss.str().c_str());
  }

  PetscFunctionReturn(0);
}

PetscErrorCode OpSetPiolaTransform::doWork(
    int side,
    EntityType type,
    DataForcesAndSurcesCore::EntData &data)  {
  PetscFunctionBegin;

  if(type != MBTRI && type != MBTET) PetscFunctionReturn(0);

  try {

  const double c = 1./6.;

  unsigned int nb_gauss_pts = data.getHdivN().size1();
  unsigned int nb_dofs = data.getHdivN().size2()/3;
  unsigned int gg = 0;
  piolaN.resize(nb_gauss_pts,data.getHdivN().size2());
  piolaDiffN.resize(nb_gauss_pts,data.getDiffHdivN().size2());
  for(;gg<nb_gauss_pts;gg++) {
    unsigned int dd = 0;
    for(;dd<nb_dofs;dd++) {
      cblas_dgemv(CblasRowMajor,CblasNoTrans,3,3,c/vOlume,
	&*Jac.data().begin(),3,&data.getHdivN()(gg,3*dd),1,0.,&piolaN(gg,3*dd),1);
      int kk = 0;
      for(;kk<3;kk++) {
	cblas_dgemv(CblasRowMajor,CblasNoTrans,3,3,c/vOlume,
	  &*Jac.data().begin(),3,&data.getDiffHdivN()(gg,9*dd+3*kk),1,0.,&piolaDiffN(gg,9*dd+3*kk),1);
      }
    }
  }
  data.getHdivN().data().swap(piolaN.data());
  data.getDiffHdivN().data().swap(piolaDiffN.data());

  } catch (exception& ex) {
    ostringstream ss;
    ss << "thorw in method: " << ex.what() << " at line " << __LINE__ << " in file " << __FILE__;
    SETERRQ(PETSC_COMM_SELF,MOFEM_STD_EXCEPTION_THROW,ss.str().c_str());
  }

  PetscFunctionReturn(0);
}


PetscErrorCode OpSetHoInvJacH1::doWork(
    int side,
    EntityType type,
    DataForcesAndSurcesCore::EntData &data)  {
  PetscFunctionBegin;

  try {

  unsigned int nb_gauss_pts = data.getN().size1();
  unsigned int nb_dofs = data.getN().size2();
  //note Vetex diffN row has size of number of gass dof
  diffNinvJac.resize(nb_gauss_pts,3*nb_dofs);
  unsigned int gg = 0;
  for(;gg<nb_gauss_pts;gg++) {
    double *inv_H = &invHoJac(gg,0);
    for(unsigned dd = 0;dd<nb_dofs;dd++) {
      double *diff_N;
      if(type == MBVERTEX) {
	diff_N = &data.getDiffN()(dd,0);
      } else {
	diff_N = &data.getDiffN()(gg,3*dd);
      }
      double *diff_N_inv_Jac = &diffNinvJac(gg,3*dd);
      cblas_dgemv(CblasRowMajor,CblasTrans,3,3,1.,inv_H,3,diff_N,1,0.,diff_N_inv_Jac,1); 
    }
  }

  if(type == MBVERTEX) {
    data.getDiffN().resize(diffNinvJac.size1(),diffNinvJac.size2());
  }
  data.getDiffN().data().swap(diffNinvJac.data());

  } catch (exception& ex) {
    ostringstream ss;
    ss << "thorw in method: " << ex.what() << " at line " << __LINE__ << " in file " << __FILE__;
    SETERRQ(PETSC_COMM_SELF,MOFEM_STD_EXCEPTION_THROW,ss.str().c_str());
  }

  PetscFunctionReturn(0);
}

PetscErrorCode OpSetHoInvJacHdiv::doWork(
    int side,
    EntityType type,
    DataForcesAndSurcesCore::EntData &data)  {
  PetscFunctionBegin;

  if(type != MBTRI && type != MBTET) PetscFunctionReturn(0);

  try {

  diffHdiv_invJac.resize(data.getDiffHdivN().size1(),data.getDiffHdivN().size2());

  unsigned int nb_gauss_pts = data.getDiffHdivN().size1();
  unsigned int nb_dofs = data.getDiffHdivN().size2()/9;

  unsigned int gg = 0;
  for(;gg<nb_gauss_pts;gg++) {
    double *inv_h = &invHoJac(gg,0);
    for(unsigned dd = 0;dd<nb_dofs;dd++) {
      const double *diff_hdiv = &(data.getDiffHdivN(gg)(dd,0));
      double *diff_hdiv_inv_jac = &diffHdiv_invJac(gg,9*dd);
      int kk = 0;
      for(;kk<3;kk++) {
	cblas_dgemv(CblasRowMajor,CblasTrans,3,3,1.,inv_h,3,&diff_hdiv[kk],3,0.,&diff_hdiv_inv_jac[kk],3); 
      }
    }
  }

  data.getDiffHdivN().data().swap(diffHdiv_invJac.data());

  } catch (exception& ex) {
    ostringstream ss;
    ss << "thorw in method: " << ex.what() << " at line " << __LINE__ << " in file " << __FILE__;
    SETERRQ(PETSC_COMM_SELF,MOFEM_STD_EXCEPTION_THROW,ss.str().c_str());
  }

  PetscFunctionReturn(0);
}


PetscErrorCode OpSetHoPiolaTransform::doWork(
  int side,EntityType type,DataForcesAndSurcesCore::EntData &data) {
  PetscFunctionBegin;

  if(type != MBTRI && type != MBTET) PetscFunctionReturn(0);

  try{

  unsigned int nb_gauss_pts = data.getHdivN().size1();
  unsigned int nb_dofs = data.getHdivN().size2()/3;
  unsigned int gg = 0;
  piolaN.resize(nb_gauss_pts,data.getHdivN().size2());
  piolaDiffN.resize(nb_gauss_pts,data.getDiffHdivN().size2());

  for(;gg<nb_gauss_pts;gg++) {
    unsigned int dd = 0;
    for(;dd<nb_dofs;dd++) {
      cblas_dgemv(CblasRowMajor,CblasNoTrans,3,3,1./detHoJac[gg],
	&hoJac(gg,0),3,&data.getHdivN()(gg,3*dd),1,0.,&piolaN(gg,3*dd),1);
      int kk = 0;
      for(;kk<3;kk++) {
	cblas_dgemv(CblasRowMajor,CblasNoTrans,3,3,1./detHoJac[gg],
	  &hoJac(gg,0),3,&data.getDiffHdivN()(gg,9*dd+3*kk),1,0.,&piolaDiffN(gg,9*dd+3*kk),1);
      }
    }
  }

  data.getHdivN().data().swap(piolaN.data());
  data.getDiffHdivN().data().swap(piolaDiffN.data());


  } catch (exception& ex) {
    ostringstream ss;
    ss << "thorw in method: " << ex.what() << " at line " << __LINE__ << " in file " << __FILE__;
    SETERRQ(PETSC_COMM_SELF,MOFEM_STD_EXCEPTION_THROW,ss.str().c_str());
  }

  PetscFunctionReturn(0);
}


PetscErrorCode OpGetData::doWork(
    int side,EntityType type,DataForcesAndSurcesCore::EntData &data) {
  PetscFunctionBegin;

  try {

  if(data.getFieldData().size() == 0) {
    PetscFunctionReturn(0);
  }

  unsigned int nb_dofs = data.getFieldData().size();
  if(nb_dofs % rank != 0) {
    SETERRQ4(PETSC_COMM_SELF,MOFEM_DATA_INSONSISTENCY,
      "data inconsistency, type %d, side %d, nb_dofs %d, rank %d",
      type,side,nb_dofs,rank);
  }
  if(nb_dofs/rank > data.getN().size2()) {
    SETERRQ2(PETSC_COMM_SELF,MOFEM_DATA_INSONSISTENCY,
      "data inconsistency nb_dofs >= data.N.size2() %u >= %u",nb_dofs,data.getN().size2());
  }
  data_at_GaussPt.resize(data.getN().size1(),rank);
  dataGrad_at_GaussPt.resize(data.getN().size1(),rank*dim);
  if(type == MBVERTEX) {
    bzero(&*data_at_GaussPt.data().begin(),data.getN().size1()*rank*sizeof(FieldData));
    bzero(&*dataGrad_at_GaussPt.data().begin(),data.getN().size1()*rank*dim*sizeof(FieldData));
    for(int rr = 0;rr<rank;rr++) {
      for(unsigned int dd = 0;dd<dim;dd++) {
	dataGrad_at_GaussPt(0,dim*rr+dd) = cblas_ddot(nb_dofs/rank,&data.getDiffN()(0,dd),dim,&data.getFieldData()[rr],rank);
      }
    }
  }
  for(unsigned int gg = 0;gg<data.getN().size1();gg++) {
    for(int rr = 0;rr<rank;rr++) {
      data_at_GaussPt(gg,rr) = cblas_ddot(nb_dofs/rank,&data.getN()(gg,0),1,&data.getFieldData()[rr],rank);
      for(unsigned int dd = 0;dd<dim;dd++) {
	if(type == MBVERTEX) {
	  if(gg == 0) continue;
	  dataGrad_at_GaussPt(gg,dim*rr+dd) += dataGrad_at_GaussPt(0,dim*rr+dd);
	} else {
	  dataGrad_at_GaussPt(gg,dim*rr+dd) += cblas_ddot(nb_dofs/rank,&data.getDiffN()(gg,dd),dim,&data.getFieldData()[rr],rank);
	}
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

PetscErrorCode TetElementForcesAndSourcesCore::operator()() {
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
    nb_gauss_pts = gm_rule_size(rule,3);
    gaussPts.resize(4,nb_gauss_pts);
    ierr = Grundmann_Moeller_integration_points_3D_TET(
      rule,&gaussPts(0,0),&gaussPts(1,0),&gaussPts(2,0),&gaussPts(3,0)); CHKERRQ(ierr);
  } else {
    ierr = setGaussPts(order); CHKERRQ(ierr);
    nb_gauss_pts = gaussPts.size2();
  }

  ierr = shapeTETFunctions_H1(dataH1,
    &gaussPts(0,0),&gaussPts(1,0),&gaussPts(2,0),nb_gauss_pts); CHKERRQ(ierr);

  if((dataH1.spacesOnEntities[MBTRI]).test(HDIV)) {
    ierr = shapeTETFunctions_Hdiv(dataHdiv,
      &gaussPts(0,0),&gaussPts(1,0),&gaussPts(2,0),nb_gauss_pts); CHKERRQ(ierr);
  }

  if((dataH1.spacesOnEntities[MBTET]).test(L2)) {
    ierr = shapeTETFunctions_L2(dataL2,
      &gaussPts(0,0),&gaussPts(1,0),&gaussPts(2,0),nb_gauss_pts); CHKERRQ(ierr);
  }

  EntityHandle ent = fePtr->get_ent();
  int num_nodes;
  const EntityHandle* conn;
  rval = mField.get_moab().get_connectivity(ent,conn,num_nodes,true); CHKERR_PETSC(rval);
  coords.resize(num_nodes*3);
  rval = mField.get_moab().get_coords(conn,num_nodes,&*coords.data().begin()); CHKERR_PETSC(rval);
  vOlume = Shape_intVolumeMBTET(&*dataH1.dataOnEntities[MBVERTEX][0].getDiffN().data().begin(),&*coords.data().begin()); 
  Jac.resize(3,3);
  invJac.resize(3,3);
  ierr = ShapeJacMBTET(&*dataH1.dataOnEntities[MBVERTEX][0].getDiffN().data().begin(),&*coords.begin(),&*Jac.data().begin()); CHKERRQ(ierr);
  noalias(invJac) = Jac;
  ierr = Shape_invJac(&*invJac.data().begin()); CHKERRQ(ierr);

  coordsAtGaussPts.resize(nb_gauss_pts,3);
  for(int gg = 0;gg<nb_gauss_pts;gg++) {
    for(int dd = 0;dd<3;dd++) {
      coordsAtGaussPts(gg,dd) = cblas_ddot(4,&dataH1.dataOnEntities[MBVERTEX][0].getN()(gg,0),1,&coords[dd],3);
    }
  }

  try {
    ierr = opSetInvJacH1.opRhs(dataH1); CHKERRQ(ierr);
    ierr = opPiolaTransform.opRhs(dataHdiv); CHKERRQ(ierr);
    ierr = opSetInvJacHdiv.opRhs(dataHdiv); CHKERRQ(ierr);
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
    try {
      ierr = opHOatGaussPoints.opRhs(dataH1); CHKERRQ(ierr);
      hoGaussPtsInvJac.resize(hoGaussPtsJac.size1(),hoGaussPtsJac.size2());
      ublas::noalias(hoGaussPtsInvJac) = hoGaussPtsJac;
      ublas::matrix<double> jac(3,3);
      hoGaussPtsDetJac.resize(nb_gauss_pts);
      for(int gg = 0;gg<nb_gauss_pts;gg++) {
	cblas_dcopy(9,&hoGaussPtsJac(gg,0),1,&jac(0,0),1);
	hoGaussPtsDetJac[gg] = Shape_detJac(&jac(0,0));
	ierr = Shape_invJac(&hoGaussPtsInvJac(gg,0)); CHKERRQ(ierr);
      }

      ierr = opSetHoInvJacH1.opRhs(dataH1); CHKERRQ(ierr);
      ierr = opSetHoPiolaTransform.opRhs(dataHdiv); CHKERRQ(ierr);
      ierr = opSetHoInvJacHdiv.opRhs(dataHdiv); CHKERRQ(ierr);

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

  for(
    boost::ptr_vector<UserDataOperator>::iterator oit = vecUserOpN.begin();
    oit != vecUserOpN.end(); oit++) {

    oit->setPtrFE(this);
    BitFieldId data_id = mField.get_field_structure(oit->row_field_name)->get_id();
    if((oit->getMoFEMFEPtr()->get_BitFieldId_data()&data_id).none()) {
      SETERRQ1(PETSC_COMM_SELF,MOFEM_DATA_INSONSISTENCY,"no data field < %s > on finite elemeny",oit->row_field_name.c_str());
    }

    FieldSpace row_space = mField.get_field_structure(oit->row_field_name)->get_space();
    
    DataForcesAndSurcesCore *op_data = NULL;
    switch(row_space) {
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
	op_data = &dataL2;
	break;
      default:
	SETERRQ(PETSC_COMM_SELF,MOFEM_DATA_INSONSISTENCY,"data inconsistency");
      break;
    }

    switch(row_space) {
      case H1:
      ierr = getRowNodesIndices(*op_data,oit->row_field_name); CHKERRQ(ierr);
      ierr = getNodesFieldData(*op_data,oit->row_field_name); CHKERRQ(ierr);
      case HCURL:
      ierr = getEdgesRowIndices(*op_data,oit->row_field_name); CHKERRQ(ierr);
      ierr = getEdgesOrder(*op_data,oit->row_field_name); CHKERRQ(ierr);
      ierr = getEdgesFieldData(*op_data,oit->row_field_name); CHKERRQ(ierr);
      case HDIV:
      ierr = getTrisRowIndices(*op_data,oit->row_field_name); CHKERRQ(ierr);
      ierr = getTrisOrder(*op_data,oit->row_field_name); CHKERRQ(ierr);
      ierr = getTrisFieldData(*op_data,oit->row_field_name); CHKERRQ(ierr);
      case L2:
      ierr = getTetsRowIndices(*op_data,oit->row_field_name); CHKERRQ(ierr);
      ierr = getTetsOrder(*op_data,oit->row_field_name); CHKERRQ(ierr);
      ierr = getTetsFieldData(*op_data,oit->row_field_name); CHKERRQ(ierr);
      default:
      break;
    }

    try {
      ierr = oit->opRhs(*op_data); CHKERRQ(ierr);
    } catch (exception& ex) {
      ostringstream ss;
      ss << "thorw in method: " << ex.what() << " at line " << __LINE__ << " in file " << __FILE__;
      SETERRQ(PETSC_COMM_SELF,MOFEM_STD_EXCEPTION_THROW,ss.str().c_str());
    }

  }

  for(
    boost::ptr_vector<UserDataOperator>::iterator oit = vecUserOpNN.begin();
    oit != vecUserOpNN.end(); oit++) {

    oit->setPtrFE(this);
    BitFieldId row_id = mField.get_field_structure(oit->row_field_name)->get_id();
    BitFieldId col_id = mField.get_field_structure(oit->col_field_name)->get_id();

    if((oit->getMoFEMFEPtr()->get_BitFieldId_data()&row_id).none()) {
      SETERRQ1(PETSC_COMM_SELF,MOFEM_NOT_FOUND,"no row field < %s > on finite elemeny",oit->row_field_name.c_str());
    }
    if((oit->getMoFEMFEPtr()->get_BitFieldId_data()&col_id).none()) {
      SETERRQ1(PETSC_COMM_SELF,MOFEM_NOT_FOUND,"no data field < %s > on finite elemeny",oit->col_field_name.c_str());
    }

    FieldSpace row_space = mField.get_field_structure(oit->row_field_name)->get_space();

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
	row_op_data = &dataL2;
	break;
      default:
	SETERRQ(PETSC_COMM_SELF,MOFEM_DATA_INSONSISTENCY,"data inconsistency");
      break;
    }

    switch(row_space) {
      case H1:
      ierr = getRowNodesIndices(*row_op_data,oit->row_field_name); CHKERRQ(ierr);
      ierr = getNodesFieldData(*row_op_data,oit->row_field_name); CHKERRQ(ierr);
      case HCURL:
      ierr = getEdgesRowIndices(*row_op_data,oit->row_field_name); CHKERRQ(ierr);
      ierr = getEdgesOrder(*row_op_data,oit->row_field_name); CHKERRQ(ierr);
      ierr = getEdgesFieldData(*row_op_data,oit->row_field_name); CHKERRQ(ierr);
      case HDIV:
      ierr = getTrisRowIndices(*row_op_data,oit->row_field_name); CHKERRQ(ierr);
      ierr = getTrisOrder(*row_op_data,oit->row_field_name); CHKERRQ(ierr);
      ierr = getTrisFieldData(*row_op_data,oit->row_field_name); CHKERRQ(ierr);
      case L2:
      ierr = getTetsRowIndices(*row_op_data,oit->row_field_name); CHKERRQ(ierr);
      ierr = getTetsOrder(*row_op_data,oit->row_field_name); CHKERRQ(ierr);
      ierr = getTetsFieldData(*row_op_data,oit->row_field_name); CHKERRQ(ierr);
      default:
      break;
    }

    FieldSpace col_space = mField.get_field_structure(oit->col_field_name)->get_space();
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
	col_op_data = &derivedDataL2;
	break;
      default:
	SETERRQ(PETSC_COMM_SELF,MOFEM_DATA_INSONSISTENCY,"data inconsistency");
      break;
    }

    switch(col_space) {
      case H1:
      ierr = getColNodesIndices(*col_op_data,oit->col_field_name); CHKERRQ(ierr);
      ierr = getNodesFieldData(*col_op_data,oit->col_field_name); CHKERRQ(ierr);
      case HCURL:
      ierr = getEdgesColIndices(*col_op_data,oit->col_field_name); CHKERRQ(ierr);
      ierr = getEdgesOrder(*col_op_data,oit->col_field_name); CHKERRQ(ierr);
      ierr = getEdgesFieldData(*col_op_data,oit->col_field_name); CHKERRQ(ierr);
      case HDIV:
      ierr = getTrisColIndices(*col_op_data,oit->col_field_name); CHKERRQ(ierr);
      ierr = getTrisOrder(*col_op_data,oit->col_field_name); CHKERRQ(ierr);
      ierr = getTrisFieldData(*col_op_data,oit->col_field_name); CHKERRQ(ierr);
      case L2:
      ierr = getTetsColIndices(*col_op_data,oit->col_field_name); CHKERRQ(ierr);
      ierr = getTetsOrder(*col_op_data,oit->col_field_name); CHKERRQ(ierr);
      ierr = getTetsFieldData(*col_op_data,oit->col_field_name); CHKERRQ(ierr);
      default:
      break;
    }

    try {
      ierr = oit->opLhs(*row_op_data,*col_op_data,oit->symm); CHKERRQ(ierr);
    } catch (exception& ex) {
      ostringstream ss;
      ss << "thorw in method: " << ex.what() << " at line " << __LINE__ << " in file " << __FILE__;
      SETERRQ(PETSC_COMM_SELF,MOFEM_STD_EXCEPTION_THROW,ss.str().c_str());
    }

  }

  } catch (exception& ex) {
    ostringstream ss;
    ss << "thorw in method: " << ex.what() << " at line " << __LINE__ << " in file " << __FILE__;
    SETERRQ(PETSC_COMM_SELF,MOFEM_STD_EXCEPTION_THROW,ss.str().c_str());
  }

  PetscFunctionReturn(0);
}

PetscErrorCode TetElementForcesAndSourcesCore::UserDataOperator::getDivergenceMatrixOperato_Hdiv(
      int side,EntityType type,DataForcesAndSurcesCore::EntData &data,
      int gg,ublas::vector<FieldData> &div) {
  PetscFunctionBegin;

  try {

  int nb_dofs = data.getFieldData().size();
  if((unsigned int)nb_dofs != data.getDiffHdivN().size2()/9) {
    SETERRQ(PETSC_COMM_SELF,MOFEM_DATA_INSONSISTENCY,"data inconsistency");
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

PetscErrorCode OpGetNormals::doWork(int side,EntityType type,DataForcesAndSurcesCore::EntData &data) {
  PetscFunctionBegin;
  
  if(data.getFieldData().size()==0)  PetscFunctionReturn(0);
  
  switch (type) {
    case MBVERTEX: {
      for(unsigned int gg = 0;gg<data.getN().size1();gg++) {
	for(int nn = 0;nn<3;nn++) {
	  tAngent1_at_GaussPt(gg,nn) = cblas_ddot(3,&data.getDiffN()(0,0),2,&data.getFieldData()[nn],3);
	  tAngent2_at_GaussPt(gg,nn) = cblas_ddot(3,&data.getDiffN()(0,1),2,&data.getFieldData()[nn],3);
	}
      }
    } 
    break;
    case MBEDGE:     
    case MBTRI: {
      /*cerr << side << " " << type << endl;
      cerr << data.getN() << endl;
      cerr << data.getDiffN() << endl;
      cerr << data.getFieldData() << endl;
      cerr << "t1 " << tAngent1_at_GaussPt << endl;
      cerr << "t2 " << tAngent2_at_GaussPt << endl;*/
      if(2*data.getN().size2() != data.getDiffN().size2()) {
	SETERRQ(PETSC_COMM_SELF,MOFEM_DATA_INSONSISTENCY,"data inconsistency");
      }
      unsigned int nb_dofs = data.getFieldData().size();
      if(nb_dofs%3!=0) {
	SETERRQ(PETSC_COMM_SELF,MOFEM_DATA_INSONSISTENCY,"data inconsistency");
      }
      if(nb_dofs > 3*data.getN().size2()) {
	SETERRQ(PETSC_COMM_SELF,MOFEM_DATA_INSONSISTENCY,"data inconsistency");
      }
      for(unsigned int gg = 0;gg<data.getN().size1();gg++) {
	for(int dd = 0;dd<3;dd++) {
	  tAngent1_at_GaussPt(gg,dd) += cblas_ddot(nb_dofs/3,&data.getDiffN()(gg,0),2,&data.getFieldData()[dd],3);
	  tAngent2_at_GaussPt(gg,dd) += cblas_ddot(nb_dofs/3,&data.getDiffN()(gg,1),2,&data.getFieldData()[dd],3);
	}
      }
      //cerr << "t1 " << tAngent1_at_GaussPt << endl;
      //cerr << "t2 " << tAngent2_at_GaussPt << endl;
    }
    break;
    default:
      SETERRQ(PETSC_COMM_SELF,MOFEM_NOT_IMPLEMENTED,"not implemented");
  }

  PetscFunctionReturn(0);
}

PetscErrorCode OpSetPiolaTransoformOnTriangle::doWork(
    int side,
    EntityType type,
    DataForcesAndSurcesCore::EntData &data) {
  PetscFunctionBegin;

  if(type != MBTRI) PetscFunctionReturn(0);
	
  double l0 = cblas_dnrm2(3,&normal[0],1);
  int nb_gauss_pts = data.getHdivN().size1();
  int nb_dofs = data.getHdivN().size2()/3;
  int gg = 0;
  for(;gg<nb_gauss_pts;gg++) {
    
    int dd = 0;
    for(;dd<nb_dofs;dd++) {
      double val = data.getHdivN()(gg,3*dd);
      if(nOrmals_at_GaussPt.size1()==(unsigned int)nb_gauss_pts) {
	double l = cblas_dnrm2(3,&nOrmals_at_GaussPt(gg,0),1);
	cblas_dcopy(3,&nOrmals_at_GaussPt(gg,0),1,&data.getHdivN()(gg,3*dd),1);
	cblas_dscal(3,val/pow(l,2),&data.getHdivN()(gg,3*dd),1);
      } else {
	cblas_dcopy(3,&normal[0],1,&data.getHdivN()(gg,3*dd),1);
	cblas_dscal(3,val/pow(l0,2),&data.getHdivN()(gg,3*dd),1);
      }
    }    

  }

  PetscFunctionReturn(0);
}

PetscErrorCode OpGetNormals::calculateNormals() {
  PetscFunctionBegin;
  PetscErrorCode ierr;

  sPin.resize(3,3);
  bzero(&*sPin.data().begin(),9*sizeof(FieldData));
  nOrmals_at_GaussPt.resize(tAngent1_at_GaussPt.size1(),3);
  for(unsigned int gg = 0;gg<tAngent1_at_GaussPt.size1();gg++) {
    ierr = Spin(&*sPin.data().begin(),&tAngent1_at_GaussPt(gg,0)); CHKERRQ(ierr);
    cblas_dgemv(
      CblasRowMajor,CblasNoTrans,3,3,1.,
      &*sPin.data().begin(),3,&tAngent2_at_GaussPt(gg,0),1,0.,
      &nOrmals_at_GaussPt(gg,0),1);
  }

  PetscFunctionReturn(0);
}

PetscErrorCode TriElementForcesAndSurcesCore::operator()() {
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
    nb_gauss_pts = gm_rule_size(rule,2);
    gaussPts.resize(3,nb_gauss_pts);
    ierr = Grundmann_Moeller_integration_points_2D_TRI(
      rule,&gaussPts(0,0),&gaussPts(1,0),&gaussPts(2,0)); CHKERRQ(ierr);
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
    &*coords.data().begin(),&*normal.data().begin()); CHKERRQ(ierr);
  aRea = cblas_dnrm2(3,&*normal.data().begin(),1)*0.5;

  coordsAtGaussPts.resize(nb_gauss_pts,3);
  for(int gg = 0;gg<nb_gauss_pts;gg++) {
    for(int dd = 0;dd<3;dd++) {
      coordsAtGaussPts(gg,dd) = cblas_ddot(3,&dataH1.dataOnEntities[MBVERTEX][0].getN()(gg,0),1,&coords[dd],3);
    }
  }

  if(mField.check_field(meshPositionsFieldName)) {
    nOrmals_at_GaussPt.resize(nb_gauss_pts,3);
    tAngent1_at_GaussPt.resize(nb_gauss_pts,3);
    tAngent2_at_GaussPt.resize(nb_gauss_pts,3);
    ierr = getEdgesOrder(dataH1,meshPositionsFieldName); CHKERRQ(ierr);
    ierr = getTrisOrder(dataH1,meshPositionsFieldName); CHKERRQ(ierr);
    ierr = getNodesFieldData(dataH1,meshPositionsFieldName); CHKERRQ(ierr);
    ierr = getEdgesFieldData(dataH1,meshPositionsFieldName); CHKERRQ(ierr);
    ierr = getTrisFieldData(dataH1,meshPositionsFieldName); CHKERRQ(ierr);
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

  for(
    boost::ptr_vector<UserDataOperator>::iterator oit = vecUserOpN.begin();
    oit != vecUserOpN.end(); oit++) {

    oit->setPtrFE(this);
    BitFieldId data_id = mField.get_field_structure(oit->row_field_name)->get_id();
    if((oit->getMoFEMFEPtr()->get_BitFieldId_data()&data_id).none()) {
      SETERRQ1(PETSC_COMM_SELF,MOFEM_DATA_INSONSISTENCY,"no data field < %s > on finite elemeny",oit->row_field_name.c_str());
    }

    FieldSpace row_space = mField.get_field_structure(oit->row_field_name)->get_space();
    
    DataForcesAndSurcesCore *op_data = NULL;
    switch(row_space) {
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
      default:
	SETERRQ(PETSC_COMM_SELF,MOFEM_DATA_INSONSISTENCY,"data inconsistency");
      break;
    }

    switch(row_space) {
      case H1:
      ierr = getRowNodesIndices(*op_data,oit->row_field_name); CHKERRQ(ierr);
      ierr = getNodesFieldData(*op_data,oit->row_field_name); CHKERRQ(ierr);
      case HCURL:
      ierr = getEdgesRowIndices(*op_data,oit->row_field_name); CHKERRQ(ierr);
      ierr = getEdgesOrder(*op_data,oit->row_field_name); CHKERRQ(ierr);
      ierr = getEdgesFieldData(*op_data,oit->row_field_name); CHKERRQ(ierr);
      case HDIV:
      case L2:
      ierr = getTrisRowIndices(*op_data,oit->row_field_name); CHKERRQ(ierr);
      ierr = getTrisOrder(*op_data,oit->row_field_name); CHKERRQ(ierr);
      ierr = getTrisFieldData(*op_data,oit->row_field_name); CHKERRQ(ierr);
      default:
      break;
    }

    try {
      ierr = oit->opRhs(*op_data); CHKERRQ(ierr);
    } catch (exception& ex) {
      ostringstream ss;
      ss << "thorw in method: " << ex.what() << " at line " << __LINE__ << " in file " << __FILE__;
      SETERRQ(PETSC_COMM_SELF,MOFEM_STD_EXCEPTION_THROW,ss.str().c_str());
    }

  }

  for(
    boost::ptr_vector<UserDataOperator>::iterator oit = vecUserOpSymmNN.begin();
    oit != vecUserOpSymmNN.end(); oit++) {

    oit->setPtrFE(this);
    BitFieldId row_id = mField.get_field_structure(oit->row_field_name)->get_id();
    BitFieldId col_id = mField.get_field_structure(oit->col_field_name)->get_id();

    if((oit->getMoFEMFEPtr()->get_BitFieldId_data()&row_id).none()) {
      SETERRQ1(PETSC_COMM_SELF,MOFEM_NOT_FOUND,"no row field < %s > on finite elemeny",oit->row_field_name.c_str());
    }
    if((oit->getMoFEMFEPtr()->get_BitFieldId_data()&col_id).none()) {
      SETERRQ1(PETSC_COMM_SELF,MOFEM_NOT_FOUND,"no data field < %s > on finite elemeny",oit->col_field_name.c_str());
    }

    FieldSpace row_space = mField.get_field_structure(oit->row_field_name)->get_space();

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
      default:
	SETERRQ(PETSC_COMM_SELF,MOFEM_DATA_INSONSISTENCY,"data inconsistency");
      break;
    }

    switch(row_space) {
      case H1:
      ierr = getRowNodesIndices(*row_op_data,oit->row_field_name); CHKERRQ(ierr);
      ierr = getNodesFieldData(*row_op_data,oit->row_field_name); CHKERRQ(ierr);
      case HCURL:
      ierr = getEdgesRowIndices(*row_op_data,oit->row_field_name); CHKERRQ(ierr);
      ierr = getEdgesOrder(*row_op_data,oit->row_field_name); CHKERRQ(ierr);
      ierr = getEdgesFieldData(*row_op_data,oit->row_field_name); CHKERRQ(ierr);
      case HDIV:
      case L2:
      ierr = getTrisRowIndices(*row_op_data,oit->row_field_name); CHKERRQ(ierr);
      ierr = getTrisOrder(*row_op_data,oit->row_field_name); CHKERRQ(ierr);
      ierr = getTrisFieldData(*row_op_data,oit->row_field_name); CHKERRQ(ierr);
      default:
      break;
    }

    FieldSpace col_space = mField.get_field_structure(oit->col_field_name)->get_space();
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
	SETERRQ(PETSC_COMM_SELF,MOFEM_DATA_INSONSISTENCY,"data inconsistency");
      break;
    }

    switch(col_space) {
      case H1:
      ierr = getColNodesIndices(*col_op_data,oit->col_field_name); CHKERRQ(ierr);
      ierr = getNodesFieldData(*col_op_data,oit->col_field_name); CHKERRQ(ierr);
      case HCURL:
      ierr = getEdgesColIndices(*col_op_data,oit->col_field_name); CHKERRQ(ierr);
      ierr = getEdgesOrder(*col_op_data,oit->col_field_name); CHKERRQ(ierr);
      ierr = getEdgesFieldData(*col_op_data,oit->col_field_name); CHKERRQ(ierr);
      case HDIV:
      case L2:
      ierr = getTrisColIndices(*col_op_data,oit->col_field_name); CHKERRQ(ierr);
      ierr = getTrisOrder(*col_op_data,oit->col_field_name); CHKERRQ(ierr);
      ierr = getTrisFieldData(*col_op_data,oit->col_field_name); CHKERRQ(ierr);
      default:
      break;
    }

    try {
      ierr = oit->opLhs(*row_op_data,*col_op_data,oit->symm); CHKERRQ(ierr);
    } catch (exception& ex) {
      ostringstream ss;
      ss << "thorw in method: " << ex.what() << " at line " << __LINE__ << " in file " << __FILE__;
      SETERRQ(PETSC_COMM_SELF,MOFEM_STD_EXCEPTION_THROW,ss.str().c_str());
    }

  }

  PetscFunctionReturn(0);
}

PetscErrorCode EdgeElementForcesAndSurcesCore::operator()() {
  PetscFunctionBegin;

  if(fePtr->get_ent_type() != MBEDGE) PetscFunctionReturn(0);

  PetscErrorCode ierr;

  //PetscAttachDebugger();

  ierr = getEdgesOrder(data,H1); CHKERRQ(ierr);

  int order = data.dataOnEntities[MBEDGE][0].getOrder();
  int rule = getRule(order);
  int nb_gauss_pts = gm_rule_size(rule,1);
  gaussPts.resize(2,nb_gauss_pts);

  ierr = Grundmann_Moeller_integration_points_1D_EDGE(rule,&gaussPts(0,0),&gaussPts(1,0)); CHKERRQ(ierr);
  ierr = shapeEDGEFunctions_H1(data,&gaussPts(0,0),nb_gauss_pts); CHKERRQ(ierr);

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

  DataForcesAndSurcesCore *col_data = &derivedData;

  for(
    boost::ptr_vector<UserDataOperator>::iterator oit = vecUserOpN.begin();
    oit != vecUserOpN.end(); oit++) {

    oit->setPtrFE(this);
    BitFieldId row_id = mField.get_field_structure(oit->row_field_name)->get_id();
    BitFieldId col_id = mField.get_field_structure(oit->col_field_name)->get_id();

    if((oit->getMoFEMFEPtr()->get_BitFieldId_row()&row_id).none()) {
      SETERRQ1(PETSC_COMM_SELF,1,"no row field < %s > on finite elemeny",oit->row_field_name.c_str());
    }
    if((oit->getMoFEMFEPtr()->get_BitFieldId_data()&col_id).none()) {
      SETERRQ1(PETSC_COMM_SELF,1,"no data field < %s > on finite elemeny",oit->col_field_name.c_str());
    }

    //row indices
    ierr = getRowNodesIndices(data,oit->row_field_name); CHKERRQ(ierr);
    ierr = getEdgesRowIndices(data,oit->row_field_name); CHKERRQ(ierr);
    //col data
    ierr = getEdgesOrder(data,oit->col_field_name); CHKERRQ(ierr);
    ierr = getNodesFieldData(data,oit->col_field_name); CHKERRQ(ierr);
    ierr = getEdgesFieldData(data,oit->col_field_name); CHKERRQ(ierr);

    try {
      ierr = oit->opRhs(data); CHKERRQ(ierr);
    } catch (exception& ex) {
      ostringstream ss;
      ss << "thorw in method: " << ex.what() << " at line " << __LINE__ << " in file " << __FILE__;
      SETERRQ(PETSC_COMM_SELF,1,ss.str().c_str());
    }

  }

  for(
    boost::ptr_vector<UserDataOperator>::iterator oit = vecUserOpSymmNN.begin();
    oit != vecUserOpSymmNN.end(); oit++) {

    oit->setPtrFE(this);
    BitFieldId row_id = mField.get_field_structure(oit->row_field_name)->get_id();
    BitFieldId col_id = mField.get_field_structure(oit->col_field_name)->get_id();

    if((oit->getMoFEMFEPtr()->get_BitFieldId_row()&row_id).none()) {
      SETERRQ1(PETSC_COMM_SELF,1,"no row field < %s > on finite elemeny",oit->row_field_name.c_str());
    }
    if((oit->getMoFEMFEPtr()->get_BitFieldId_col()&col_id).none()) {
      SETERRQ1(PETSC_COMM_SELF,1,"no data field < %s > on finite elemeny",oit->col_field_name.c_str());
    }
    if((oit->getMoFEMFEPtr()->get_BitFieldId_data()&col_id).none()) {
      SETERRQ1(PETSC_COMM_SELF,1,"no data field < %s > on finite elemeny",oit->col_field_name.c_str());
    }

    //row indices
    ierr = getEdgesOrder(data,oit->row_field_name); CHKERRQ(ierr);
    ierr = getRowNodesIndices(data,oit->row_field_name); CHKERRQ(ierr);
    ierr = getEdgesRowIndices(data,oit->row_field_name); CHKERRQ(ierr);
    //col indices
    ierr = getEdgesOrder(*col_data,oit->col_field_name); CHKERRQ(ierr);
    ierr = getColNodesIndices(*col_data,oit->col_field_name); CHKERRQ(ierr);
    //col data
    ierr = getEdgesColIndices(*col_data,oit->col_field_name); CHKERRQ(ierr);
    ierr = getNodesFieldData(data,oit->col_field_name); CHKERRQ(ierr);
    ierr = getEdgesFieldData(data,oit->col_field_name); CHKERRQ(ierr);

    try {
      ierr = oit->opLhs(data,*col_data,true); CHKERRQ(ierr);
    } catch (exception& ex) {
      ostringstream ss;
      ss << "thorw in method: " << ex.what() << " at line " << __LINE__ << " in file " << __FILE__;
      SETERRQ(PETSC_COMM_SELF,1,ss.str().c_str());
    }

  }

  PetscFunctionReturn(0);
}

PetscErrorCode VertexElementForcesAndSurcesCore::operator()() {
  PetscFunctionBegin;

  if(fePtr->get_ent_type() != MBVERTEX) PetscFunctionReturn(0);

  EntityHandle ent = fePtr->get_ent();
  coords.resize(3);
  rval = mField.get_moab().get_coords(&ent,1,&*coords.data().begin()); CHKERR_PETSC(rval);

  DataForcesAndSurcesCore *col_data = &derivedData;

  for(
    boost::ptr_vector<UserDataOperator>::iterator oit = vecUserOpN.begin();
    oit != vecUserOpN.end(); oit++) {

    oit->setPtrFE(this);
    BitFieldId row_id = mField.get_field_structure(oit->row_field_name)->get_id();
    BitFieldId col_id = mField.get_field_structure(oit->col_field_name)->get_id();

    if((oit->getMoFEMFEPtr()->get_BitFieldId_row()&row_id).none()) {
      SETERRQ1(PETSC_COMM_SELF,1,"no row field < %s > on finite elemeny",oit->row_field_name.c_str());
    }
    if((oit->getMoFEMFEPtr()->get_BitFieldId_data()&col_id).none()) {
      SETERRQ1(PETSC_COMM_SELF,1,"no data field < %s > on finite elemeny",oit->row_field_name.c_str());
    }

    ierr = getRowNodesIndices(data,oit->row_field_name); CHKERRQ(ierr);
    ierr = getNodesFieldData(data,oit->col_field_name); CHKERRQ(ierr);

    try {
      ierr = oit->opRhs(data); CHKERRQ(ierr);
    } catch (exception& ex) {
      ostringstream ss;
      ss << "thorw in method: " << ex.what() << " at line " << __LINE__ << " in file " << __FILE__;
      SETERRQ(PETSC_COMM_SELF,1,ss.str().c_str());
    }

  }

  for(
    boost::ptr_vector<UserDataOperator>::iterator oit = vecUserOpSymmNN.begin();
    oit != vecUserOpSymmNN.end(); oit++) {

    oit->setPtrFE(this);
    BitFieldId row_id = mField.get_field_structure(oit->row_field_name)->get_id();
    BitFieldId col_id = mField.get_field_structure(oit->col_field_name)->get_id();

    if((oit->getMoFEMFEPtr()->get_BitFieldId_row()&row_id).none()) {
      SETERRQ1(PETSC_COMM_SELF,1,"no row field < %s > on finite elemeny",oit->row_field_name.c_str());
    }
    if((oit->getMoFEMFEPtr()->get_BitFieldId_col()&col_id).none()) {
      SETERRQ1(PETSC_COMM_SELF,1,"no data field < %s > on finite elemeny",oit->col_field_name.c_str());
    }
    if((oit->getMoFEMFEPtr()->get_BitFieldId_data()&col_id).none()) {
      SETERRQ1(PETSC_COMM_SELF,1,"no data field < %s > on finite elemeny",oit->col_field_name.c_str());
    }

    ierr = getRowNodesIndices(data,oit->row_field_name); CHKERRQ(ierr);
    ierr = getColNodesIndices(*col_data,oit->col_field_name); CHKERRQ(ierr);
    ierr = getNodesFieldData(*col_data,oit->col_field_name); CHKERRQ(ierr);

    try {
      ierr = oit->opLhs(data,*col_data,true); CHKERRQ(ierr);
    } catch (exception& ex) {
      ostringstream ss;
      ss << "thorw in method: " << ex.what() << " at line " << __LINE__ << " in file " << __FILE__;
      SETERRQ(PETSC_COMM_SELF,1,ss.str().c_str());
    }

  }

  PetscFunctionReturn(0);
}



}
