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
extern "C" {
#include "gm_rule.h"
}

using namespace boost::numeric;

namespace MoFEM {

template<class T> 
void cOnstructor(DataForcesAndSurcesCore *data,EntityType type,T) {
    
    switch (type) {
      case MBTET:
	data->nOdes.push_back(new T());
	for(int ee = 0;ee<6;ee++) {
	  data->eDges.push_back(new T());
	}
	for(int ff = 0;ff<4;ff++) {
	  data->fAces.push_back(new T());
	}
	data->vOlumes.push_back(new T());
      break;
      default:
	throw("not implemenyed");
    }

}

DataForcesAndSurcesCore::DataForcesAndSurcesCore(EntityType type) {
  cOnstructor(this,type,EntData());
}


DerivedDataForcesAndSurcesCore::DerivedDataForcesAndSurcesCore(DataForcesAndSurcesCore &data): DataForcesAndSurcesCore() {

    boost::ptr_vector<EntData>::iterator it;
    for(it = data.nOdes.begin();it!=data.nOdes.end();it++) {
      nOdes.push_back(new DerivedEntData(*it));
    }
    for(it = data.eDges.begin();it!=data.eDges.end();it++) {
      eDges.push_back(new DerivedEntData(*it));
    }
    for(it = data.fAces.begin();it!=data.fAces.end();it++) {
      fAces.push_back(new DerivedEntData(*it));
    }
    for(it = data.vOlumes.begin();it!=data.vOlumes.end();it++) {
      vOlumes.push_back(new DerivedEntData(*it));
    }
  }

ostream& operator<<(ostream& os,const DataForcesAndSurcesCore::EntData &e) {
  os << 
    "sEnse: " << e.getSense() << endl << 
    "oRder: " << e.getOrder() << endl <<
    "iNdices: " << e.getIndices() << endl;
  os.precision(2);
  os << 
    "fieldData: " << e.getFieldData() << endl;
  os <<
    "N: " << e.getN() << endl <<
    "diffN: " << e.getDiffN();
  return os;
}

ostream& operator<<(ostream& os,const DataForcesAndSurcesCore &e) {
  for(unsigned int nn = 0;nn < e.nOdes.size(); nn++) {
    os << "nOdes[" << nn << "]" << endl << e.nOdes[nn] << endl;
  }
  for(unsigned int ee = 0;ee < e.eDges.size(); ee++) {
    os << "eDges[" << ee << "]" << endl << e.eDges[ee] << endl;
  }
  for(unsigned int ff = 0;ff < e.fAces.size(); ff++) {
    os << "fAces[" << ff << "] " << endl << e.fAces[ff] << endl;
  }
  for(unsigned int vv = 0;vv < e.vOlumes.size(); vv++) {
    os << "vOlumes[" << vv << "]" << endl << e.vOlumes[vv] << endl;
  }
  return os;
}

PetscErrorCode ForcesAndSurcesCore::getSense(EntityType type,boost::ptr_vector<DataForcesAndSurcesCore::EntData> &data) {
    PetscFunctionBegin;
    SideNumber_multiIndex& side_table = const_cast<SideNumber_multiIndex&>(fe_ptr->get_side_number_table());
    if(data.size() != side_table.get<2>().count(type)) {
      SETERRQ2(PETSC_COMM_SELF,1,"data inconsistency %u != %u",data.size(),side_table.get<2>().count(type));
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
    ierr = getSense(MBEDGE,data.eDges); CHKERRQ(ierr);
    PetscFunctionReturn(0);
  }

PetscErrorCode ForcesAndSurcesCore::getFacesSense(DataForcesAndSurcesCore &data) {
    PetscFunctionBegin;
    PetscErrorCode ierr;
    ierr = getSense(MBTRI,data.fAces); CHKERRQ(ierr);
    PetscFunctionReturn(0);
  }

PetscErrorCode ForcesAndSurcesCore::getOrder(EntityType type,boost::ptr_vector<DataForcesAndSurcesCore::EntData> &data) {
    PetscFunctionBegin;
    SideNumber_multiIndex& side_table = const_cast<SideNumber_multiIndex&>(fe_ptr->get_side_number_table());
    if(data.size() != side_table.get<2>().count(type)) {
      SETERRQ(PETSC_COMM_SELF,1,"data inconsistency");
    }
    FEDofMoFEMEntity_multiIndex::index<EntType_mi_tag>::type &data_dofs =
      const_cast<FEDofMoFEMEntity_multiIndex::index<EntType_mi_tag>::type&>(fe_ptr->get_data_dofs().get<EntType_mi_tag>());
    FEDofMoFEMEntity_multiIndex::index<EntType_mi_tag>::type::iterator dit,hi_dit;
    dit = data_dofs.lower_bound(type);
    hi_dit = data_dofs.upper_bound(type);
    for(;dit!=hi_dit;dit++) {
      ApproximationOrder ent_order = dit->get_max_order();
      int side_number = dit->side_number_ptr->side_number;
      data[side_number].getOrder() = data[side_number].getOrder() > ent_order ? data[side_number].getOrder() : ent_order;
    }
    PetscFunctionReturn(0);
  }


PetscErrorCode ForcesAndSurcesCore::getEdgesOrder(DataForcesAndSurcesCore &data) {
    PetscFunctionBegin;
    PetscErrorCode ierr;
    ierr = getOrder(MBEDGE,data.eDges); CHKERRQ(ierr);
    PetscFunctionReturn(0);
  }

PetscErrorCode ForcesAndSurcesCore::getFacesOrder(DataForcesAndSurcesCore &data) {
    PetscFunctionBegin;
    PetscErrorCode ierr;
    ierr = getOrder(MBTRI,data.fAces); CHKERRQ(ierr);
    PetscFunctionReturn(0);
  }

PetscErrorCode ForcesAndSurcesCore::getOrderVolume(DataForcesAndSurcesCore &data) {
    PetscFunctionBegin;
    FEDofMoFEMEntity_multiIndex::index<EntType_mi_tag>::type &data_dofs =
      const_cast<FEDofMoFEMEntity_multiIndex::index<EntType_mi_tag>::type&>(fe_ptr->get_data_dofs().get<EntType_mi_tag>());
    FEDofMoFEMEntity_multiIndex::index<EntType_mi_tag>::type::iterator dit,hi_dit;
    dit = data_dofs.lower_bound(MBTET);
    hi_dit = data_dofs.upper_bound(MBTET);
    data.vOlumes[0].getOrder() = -1;
    for(;dit!=hi_dit;dit++) {
      ApproximationOrder ent_order = dit->get_max_order();
      int side_number = dit->side_number_ptr->side_number;
      data.vOlumes[0].getOrder() = data.vOlumes[0].getOrder() > ent_order ? data.vOlumes[0].getOrder() : ent_order;
    }
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
      int field_rank = dit->get_max_rank();
      nodes_indices[side_number*dit->get_max_rank()+dit->get_dof_rank()] = idx;
    }
    PetscFunctionReturn(0);
  }

PetscErrorCode ForcesAndSurcesCore::getRowNodesIndices(DataForcesAndSurcesCore &data,const string &field_name) {
    PetscFunctionBegin;
    PetscErrorCode ierr;
    ierr = getNodesIndices(
      field_name,const_cast<FENumeredDofMoFEMEntity_multiIndex&>(fe_ptr->get_rows_dofs()),data.nOdes[0].getIndices()); CHKERRQ(ierr);
    PetscFunctionReturn(0);
  }

PetscErrorCode ForcesAndSurcesCore::getColNodesIndices(DataForcesAndSurcesCore &data,const string &field_name) {
    PetscFunctionBegin;
    PetscErrorCode ierr;
    ierr = getNodesIndices(field_name,const_cast<FENumeredDofMoFEMEntity_multiIndex&>(fe_ptr->get_cols_dofs()),data.nOdes[0].getIndices()); CHKERRQ(ierr);
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
      int side_number = dit->side_number_ptr->side_number;
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
    SideNumber_multiIndex& side_table = const_cast<SideNumber_multiIndex&>(fe_ptr->get_side_number_table());
    SideNumber_multiIndex::nth_index<2>::type::iterator siit = side_table.get<2>().lower_bound(type);
    SideNumber_multiIndex::nth_index<2>::type::iterator hi_siit = side_table.get<2>().upper_bound(type);
    for(;siit!=hi_siit;siit++) {
      ierr = getTypeIndices(field_name,dofs,type,siit->side_number,data[siit->side_number].getIndices()); CHKERRQ(ierr);
    }
    PetscFunctionReturn(0);
  }

PetscErrorCode ForcesAndSurcesCore::getEdgeRowIndices(DataForcesAndSurcesCore &data,const string &field_name) {
    PetscFunctionBegin;
    PetscErrorCode ierr;
    ierr = getTypeIndices(
      field_name,const_cast<FENumeredDofMoFEMEntity_multiIndex&>(fe_ptr->get_rows_dofs()),MBEDGE,data.eDges); CHKERRQ(ierr);
    PetscFunctionReturn(0);
  }

PetscErrorCode ForcesAndSurcesCore::getEdgeColIndices(DataForcesAndSurcesCore &data,const string &field_name) {
    PetscFunctionBegin;
    PetscErrorCode ierr;
    ierr = getTypeIndices(
      field_name,const_cast<FENumeredDofMoFEMEntity_multiIndex&>(fe_ptr->get_cols_dofs()),MBEDGE,data.eDges); CHKERRQ(ierr);
    PetscFunctionReturn(0);
  }

PetscErrorCode ForcesAndSurcesCore::getFacesRowIndices(
    DataForcesAndSurcesCore &data,const string &field_name) {
    PetscFunctionBegin;
    PetscErrorCode ierr;
    ierr = getTypeIndices(field_name,const_cast<FENumeredDofMoFEMEntity_multiIndex&>(fe_ptr->get_rows_dofs()),MBTRI,data.fAces); CHKERRQ(ierr);
    PetscFunctionReturn(0);
  }

PetscErrorCode ForcesAndSurcesCore::getFacesColIndices(DataForcesAndSurcesCore &data,const string &field_name) {
    PetscFunctionBegin;
    PetscErrorCode ierr;
    ierr = getTypeIndices(field_name,const_cast<FENumeredDofMoFEMEntity_multiIndex&>(fe_ptr->get_cols_dofs()),MBTRI,data.fAces); CHKERRQ(ierr);
    PetscFunctionReturn(0);
  }

PetscErrorCode ForcesAndSurcesCore::getTetRowIndices(DataForcesAndSurcesCore &data,const string &field_name) {
    PetscFunctionBegin;
    PetscErrorCode ierr;
    ierr = getTypeIndices(field_name,const_cast<FENumeredDofMoFEMEntity_multiIndex&>(fe_ptr->get_rows_dofs()),MBTET,0,data.vOlumes[0].getIndices()); CHKERRQ(ierr);
    PetscFunctionReturn(0);
  }

PetscErrorCode ForcesAndSurcesCore::getTetColIndices(DataForcesAndSurcesCore &data,const string &field_name) {
    PetscFunctionBegin;
    PetscErrorCode ierr;
    ierr = getTypeIndices(field_name,const_cast<FENumeredDofMoFEMEntity_multiIndex&>(fe_ptr->get_cols_dofs()),MBTET,0,data.vOlumes[0].getIndices()); CHKERRQ(ierr);
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
      int field_rank = dit->get_max_rank();
      nodes_field_data[side_number*dit->get_max_rank()+dit->get_dof_rank()] = val;
    }
    PetscFunctionReturn(0);
  }

PetscErrorCode ForcesAndSurcesCore::getNodesFieldData(DataForcesAndSurcesCore &data,const string &field_name) {
    PetscFunctionBegin;
    PetscErrorCode ierr;
    ierr = getNodesFieldData(
      field_name,const_cast<FEDofMoFEMEntity_multiIndex&>(fe_ptr->get_data_dofs()),data.nOdes[0].getFieldData()); CHKERRQ(ierr);
    PetscFunctionReturn(0);
  }

PetscErrorCode ForcesAndSurcesCore::getTypeFieldData(
  const string &field_name,FEDofMoFEMEntity_multiIndex &dofs,EntityType type,int side_number,ublas::vector<FieldData> &ent_field_data) {
    PetscFunctionBegin;
    FEDofMoFEMEntity_multiIndex::index<Composite_Name_Type_And_Side_Number_mi_tag>::type::iterator dit,hi_dit;
    dit = dofs.get<Composite_Name_Type_And_Side_Number_mi_tag>().lower_bound(boost::make_tuple(field_name,type,side_number));
    hi_dit = dofs.get<Composite_Name_Type_And_Side_Number_mi_tag>().upper_bound(boost::make_tuple(field_name,type,side_number));
    for(;dit!=hi_dit;dit++) {
      FieldData val = dit->get_FieldData();
      int side_number = dit->side_number_ptr->side_number;
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
    SideNumber_multiIndex& side_table = const_cast<SideNumber_multiIndex&>(fe_ptr->get_side_number_table());
    if(data.size() != side_table.get<2>().count(type)) {
      SETERRQ(PETSC_COMM_SELF,1,"data inconsistency");
    }
    SideNumber_multiIndex::nth_index<2>::type::iterator siit = side_table.get<2>().lower_bound(type);
    SideNumber_multiIndex::nth_index<2>::type::iterator hi_siit = side_table.get<2>().upper_bound(type);
    for(;siit!=hi_siit;siit++) {
      ierr = getTypeFieldData(field_name,dofs,type,siit->side_number,data[siit->side_number].getFieldData()); CHKERRQ(ierr);
    }
    PetscFunctionReturn(0);
  }

PetscErrorCode ForcesAndSurcesCore::getEdgeFieldData(DataForcesAndSurcesCore &data,const string &field_name) {
    PetscFunctionBegin;
    PetscErrorCode ierr;
    ierr = getTypeFieldData(
      field_name,const_cast<FEDofMoFEMEntity_multiIndex&>(fe_ptr->get_data_dofs()),MBEDGE,data.eDges); CHKERRQ(ierr);
    PetscFunctionReturn(0);
  }

PetscErrorCode ForcesAndSurcesCore::getFacesFieldData(
    DataForcesAndSurcesCore &data,const string &field_name) {
    PetscFunctionBegin;
    PetscErrorCode ierr;
    ierr = getTypeFieldData(field_name,const_cast<FEDofMoFEMEntity_multiIndex&>(fe_ptr->get_data_dofs()),MBTRI,data.fAces); CHKERRQ(ierr);
    PetscFunctionReturn(0);
  }

PetscErrorCode ForcesAndSurcesCore::getTetFieldData(DataForcesAndSurcesCore &data,const string &field_name) {
    PetscFunctionBegin;
    PetscErrorCode ierr;
    ierr = getTypeFieldData(field_name,const_cast<FEDofMoFEMEntity_multiIndex&>(fe_ptr->get_data_dofs()),MBTET,0,data.vOlumes[0].getFieldData()); CHKERRQ(ierr);
    PetscFunctionReturn(0);
  }

PetscErrorCode ForcesAndSurcesCore::getFaceNodes(DataForcesAndSurcesCore &data) {
    PetscFunctionBegin;
    //PetscAttachDebugger();
    ErrorCode rval;
    data.facesNodes.resize(4,3);
    SideNumber_multiIndex& side_table = const_cast<SideNumber_multiIndex&>(fe_ptr->get_side_number_table());
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
	EntityHandle ent = fe_ptr->get_ent();
	rval = mField.get_moab().get_connectivity(ent,conn_tet,num_nodes_tet,true); CHKERR_PETSC(rval);
	if(num_nodes_tet != 4) SETERRQ(PETSC_COMM_SELF,1,"data inconsistency");
	int num_nodes_face;
	const EntityHandle *conn_face;
	rval = mField.get_moab().get_connectivity(side->ent,conn_face,num_nodes_face,true); CHKERR_PETSC(rval);
	if(num_nodes_face != 3) SETERRQ(PETSC_COMM_SELF,1,"data inconsistency");
	if(conn_face[0] != conn_tet[face_conn[0]]) SETERRQ(PETSC_COMM_SELF,1,"data inconsistency");
	if(conn_face[1] != conn_tet[face_conn[1]]) SETERRQ(PETSC_COMM_SELF,1,"data inconsistency");
	if(conn_face[2] != conn_tet[face_conn[2]]) SETERRQ(PETSC_COMM_SELF,1,"data inconsistency");
      }
    }
    PetscFunctionReturn(0);
  }

PetscErrorCode ForcesAndSurcesCore::shapeTETFunctions_H1(
    DataForcesAndSurcesCore &data,
    const double *G_X,const double *G_Y,const double *G_Z,const int G_DIM) {
    PetscFunctionBegin;

    PetscErrorCode ierr;

    data.nOdes[0].getN().resize(G_DIM,4);
    ierr = ShapeMBTET(&*data.nOdes[0].getN().data().begin(),G_X,G_Y,G_Z,G_DIM); CHKERRQ(ierr);
    data.nOdes[0].getDiffN().resize(4,3);
    ierr = ShapeDiffMBTET(&*data.nOdes[0].getDiffN().data().begin()); CHKERRQ(ierr);

    //edges
    if(data.eDges.size()!=6) {
      SETERRQ(PETSC_COMM_SELF,1,"data inconsistency");
    }
    int _sense_[6],_order_[6];
    double *_H1edgeN_[6],*_diffH1edgeN_[6];
    for(int ee = 0;ee<6;ee++) {
      if(data.eDges[ee].getSense() == 0) {
	SETERRQ(PETSC_COMM_SELF,1,"data inconsistency");
      }
      _sense_[ee] = data.eDges[ee].getSense();
      _order_[ee] = data.eDges[ee].getOrder();
      int nb_dofs = NBEDGE_H1(data.eDges[ee].getOrder());
      data.eDges[ee].getN().resize(G_DIM,nb_dofs);
      data.eDges[ee].getDiffN().resize(G_DIM,3*nb_dofs);
      _H1edgeN_[ee] = &*data.eDges[ee].getN().data().begin();
      _diffH1edgeN_[ee] = &*data.eDges[ee].getDiffN().data().begin();
    }
    ierr = H1_EdgeShapeFunctions_MBTET(
      _sense_,_order_,
      &*data.nOdes[0].getN().data().begin(),data.nOdes[0].getDiffN().data().begin(),
      _H1edgeN_,_diffH1edgeN_,G_DIM); CHKERRQ(ierr);

    //faces
    if(data.fAces.size()!=4) {
      SETERRQ(PETSC_COMM_SELF,1,"data inconsistency");
    }
    double *_H1faceN_[4],*_diffH1faceN_[4];
    for(int ff = 0;ff<4;ff++) {
      if(data.fAces[ff].getSense() == 0) {
	SETERRQ(PETSC_COMM_SELF,1,"data inconsistency");
      }
      int nb_dofs = NBFACE_H1(data.fAces[ff].getOrder());
      _order_[ff] = data.fAces[ff].getOrder();
      data.fAces[ff].getN().resize(G_DIM,nb_dofs);
      data.fAces[ff].getDiffN().resize(G_DIM,3*nb_dofs);
      _H1faceN_[ff] = &*data.fAces[ff].getN().data().begin();
      _diffH1faceN_[ff] = &*data.fAces[ff].getDiffN().data().begin();
    }
    if(data.facesNodes.size1() != 4) {
      SETERRQ(PETSC_COMM_SELF,1,"data inconsistency");
    }
    if(data.facesNodes.size2() != 3) {
      SETERRQ(PETSC_COMM_SELF,1,"data inconsistency");
    }
    ierr = H1_FaceShapeFunctions_MBTET(
      &*data.facesNodes.data().begin(),_order_,
      &*data.nOdes[0].getN().data().begin(),data.nOdes[0].getDiffN().data().begin(),
      _H1faceN_,_diffH1faceN_,G_DIM); CHKERRQ(ierr);

    //volume
    int nb_vol_dofs = NBVOLUME_H1(data.vOlumes[0].getOrder());
    data.vOlumes[0].getN().resize(G_DIM,nb_vol_dofs);
    data.vOlumes[0].getDiffN().resize(G_DIM,3*nb_vol_dofs);
    ierr = H1_VolumeShapeFunctions_MBTET(
      data.vOlumes[0].getOrder(),&*data.nOdes[0].getN().data().begin(),data.nOdes[0].getDiffN().data().begin(),
      &*data.vOlumes[0].getN().data().begin(),&*data.vOlumes[0].getDiffN().data().begin(),G_DIM); CHKERRQ(ierr);

    PetscFunctionReturn(0);
  }
  
PetscErrorCode ForcesAndSurcesCore::shapeTETFunctions_L2(
    DataForcesAndSurcesCore &data,
    const double *G_X,const double *G_Y,const double *G_Z,const int G_DIM) {
    PetscFunctionBegin;

    PetscErrorCode ierr;

    data.nOdes[0].getN().resize(G_DIM,4);
    ierr = ShapeMBTET(&*data.nOdes[0].getN().data().begin(),G_X,G_Y,G_Z,G_DIM); CHKERRQ(ierr);
    data.nOdes[0].getDiffN().resize(4,3);
    ierr = ShapeDiffMBTET(&*data.nOdes[0].getDiffN().data().begin()); CHKERRQ(ierr);

    data.vOlumes[0].getN().resize(G_DIM,NBVOLUME_L2(data.vOlumes[0].getOrder()));
    ierr = L2_ShapeFunctions_MBTET(
      data.vOlumes[0].getOrder(),
      &*data.nOdes[0].getN().data().begin(),&*data.nOdes[0].getDiffN().data().begin(),
      &*data.vOlumes[0].getN().data().begin(),NULL,G_DIM); CHKERRQ(ierr);
    PetscFunctionReturn(0);
  }

PetscErrorCode ForcesAndSurcesCore::shapeTRIFunctions_H1(
    DataForcesAndSurcesCore &data,
    ApproximationOrder order,const double *G_X,const double *G_Y,const int G_DIM) {
    PetscFunctionBegin;
    
    SETERRQ(PETSC_COMM_SELF,1,"not implemented");

    PetscFunctionReturn(0);
  }


PetscErrorCode DataOperator::opNH1NH1(
    DataForcesAndSurcesCore &row_data,DataForcesAndSurcesCore &col_data) {
  PetscFunctionBegin;
  PetscErrorCode ierr;

  int G_DIM = row_data.nOdes[0].getN().size1();
  ierr = doWork(
    -1,-1,MBVERTEX,MBVERTEX,
    row_data.nOdes[0],col_data.nOdes[0]); CHKERRQ(ierr);

  for(int ee = 0;ee<row_data.eDges.size();ee++) {
    int G_DIM = row_data.eDges[ee].getN().size1();
    int G_DIM_NODES = col_data.nOdes[0].getN().size1();
    if(G_DIM != G_DIM_NODES) {
      SETERRQ(PETSC_COMM_SELF,1,"data inconsistency");
    }
    ierr = doWork(
	ee,-1,MBEDGE,MBVERTEX,
	row_data.eDges[ee],col_data.nOdes[0]); CHKERRQ(ierr);
    int G_DIM_VOLUME = col_data.vOlumes[0].getN().size1();
    if(G_DIM != G_DIM_VOLUME) {
      SETERRQ(PETSC_COMM_SELF,1,"data inconsistency");
    }
    ierr = doWork(
	ee,-1,MBEDGE,MBTET,
	row_data.eDges[ee],col_data.vOlumes[0]); CHKERRQ(ierr);
    for(int EE = 0;EE<col_data.eDges.size();EE++) {
      int G_DIM_EE = col_data.eDges[EE].getN().size1();
      if(G_DIM != G_DIM_EE) {
	SETERRQ(PETSC_COMM_SELF,1,"data inconsistency");
      }
      ierr = doWork(
	  ee,EE,MBEDGE,MBEDGE,
	  row_data.eDges[ee],col_data.eDges[EE]); CHKERRQ(ierr);
    }
    for(int FF = 0;FF<col_data.fAces.size();FF++) {
      int G_DIM_FF = col_data.fAces[FF].getN().size1();
      if(G_DIM != G_DIM_FF) {
	SETERRQ(PETSC_COMM_SELF,1,"data inconsistency");
      }
      ierr = doWork(
	  ee,FF,MBEDGE,MBTRI,
	  row_data.eDges[ee],col_data.fAces[FF]); CHKERRQ(ierr);
    }
  }

  for(int ff = 0;ff<row_data.fAces.size();ff++) {
    int G_DIM = row_data.fAces[ff].getN().size1();
    int G_DIM_NODES = col_data.nOdes[0].getN().size1();
    if(G_DIM != G_DIM_NODES) {
      SETERRQ(PETSC_COMM_SELF,1,"data inconsistency");
    }
    ierr = doWork(
	ff,-1,MBTRI,MBVERTEX,
	row_data.fAces[ff],col_data.nOdes[0]); CHKERRQ(ierr);
    int G_DIM_VOLUME = col_data.vOlumes[0].getN().size1();
    if(G_DIM != G_DIM_VOLUME) {
      SETERRQ(PETSC_COMM_SELF,1,"data inconsistency");
    }
    ierr = doWork(
	ff,-1,MBTRI,MBTET,
	row_data.fAces[ff],col_data.vOlumes[0]); CHKERRQ(ierr);
    for(int EE = 0;EE<col_data.eDges.size();EE++) {
      int G_DIM_EE = col_data.eDges[EE].getN().size1();
      if(G_DIM != G_DIM_EE) {
	SETERRQ(PETSC_COMM_SELF,1,"data inconsistency");
      }
      ierr = doWork(
	  ff,EE,MBTRI,MBEDGE,
	  row_data.fAces[ff],col_data.eDges[EE]); CHKERRQ(ierr);
    }
    for(int FF = 0;FF<col_data.fAces.size();FF++) {
      int G_DIM_FF = col_data.fAces[FF].getN().size1();
      if(G_DIM != G_DIM_FF) {
	SETERRQ(PETSC_COMM_SELF,1,"data inconsistency");
      }
      ierr = doWork(
	  ff,FF,MBTRI,MBTRI,
	  row_data.fAces[ff],col_data.fAces[FF]); CHKERRQ(ierr);
    }
  }


  {

    int G_DIM = row_data.vOlumes[0].getN().size1();
    ierr = doWork(
	-1,-1,MBTET,MBTET,
	row_data.vOlumes[0],col_data.vOlumes[0]); CHKERRQ(ierr);

  }

  PetscFunctionReturn(0);
}

PetscErrorCode DataOperator::opNH1(DataForcesAndSurcesCore &data) {
  PetscFunctionBegin;
  PetscErrorCode ierr;

  ierr = doWork(-1,MBVERTEX,data.nOdes[0]); CHKERRQ(ierr);
  for(int ee = 0;ee<data.eDges.size();ee++) {
    ierr = doWork(ee,MBEDGE,data.eDges[ee]); CHKERRQ(ierr);
  }
  for(int ff = 0;ff<data.fAces.size();ff++) {
    ierr = doWork(ff,MBTRI,data.fAces[ff]); CHKERRQ(ierr);
  }
  ierr = doWork(-1,MBTET,data.vOlumes[0]); CHKERRQ(ierr);

  PetscFunctionReturn(0);
}

PetscErrorCode OpSetJac::doWork(
    int side,
    EntityType type,
    DataForcesAndSurcesCore::EntData &data) {
  PetscFunctionBegin;
  PetscErrorCode ierr;

  diffNinvJac.resize(data.getDiffN().size1(),data.getDiffN().size2());
  int nb_gauss_pts = data.getN().size1();
  int nb_dofs = data.getN().size2();
  if(type!=MBVERTEX) {
    if(nb_dofs != data.getDiffN().size2()/3) {
      SETERRQ2(PETSC_COMM_SELF,1,
	"data inconsistency nb_dofs != data.diffN.size2()/3 ( %u != %u/3 )",
	nb_dofs,data.getDiffN().size2());
    }
  }

  switch (type) {

    case MBVERTEX: {
      ierr = ShapeDiffMBTETinvJ(
	&*data.getDiffN().data().begin(),&*invJac.data().begin(),&*diffNinvJac.data().begin()); CHKERRQ(ierr);
      data.getDiffN().data().swap(diffNinvJac.data());
    }
    break;
    case MBEDGE:
    case MBTRI:
    case MBTET: {
      for(int gg = 0;gg<nb_gauss_pts;gg++) {
	for(int dd = 0;dd<nb_dofs;dd++) {
	  cblas_dgemv(CblasRowMajor,CblasTrans,3,3,1.,
	    &*invJac.data().begin(),3,&data.getDiffN()(gg,3*dd),1,0.,&diffNinvJac(gg,3*dd),1); 
	}
      }
      data.getDiffN().data().swap(diffNinvJac.data());
    }
    break;
    default:
      SETERRQ(PETSC_COMM_SELF,1,"data inconsistency");

  }

  PetscFunctionReturn(0);
}

PetscErrorCode OpGetData::doWork(
    int side,EntityType type,DataForcesAndSurcesCore::EntData &data) {
  PetscFunctionBegin;
  PetscErrorCode ierr;

  if(data.getFieldData().size() == 0) {
    SETERRQ(PETSC_COMM_SELF,1,"no data");
  }

  int nb_dofs = data.getFieldData().size();
  if(nb_dofs > data.getN().size2()) {
    SETERRQ2(PETSC_COMM_SELF,1,
      "data inconsistency nb_dofs >= data.N.size2() %u >= %u",nb_dofs,data.getN().size2());
  }
  if(nb_dofs % rank != 0) {
    SETERRQ(PETSC_COMM_SELF,1,"data inconsistency");
  }
  data_at_GaussPt.resize(data.getN().size1(),rank);
  dataGrad_at_GaussPt.resize(data.getN().size1(),rank*dim);
  if(type == MBVERTEX) {
    bzero(&*data_at_GaussPt.data().begin(),data.getN().size1()*rank*sizeof(FieldData));
    bzero(&*dataGrad_at_GaussPt.data().begin(),data.getN().size1()*rank*dim*sizeof(FieldData));
  }
  for(int gg = 0;gg<data.getN().size1();gg++) {
    for(int rr = 0;rr<rank;rr++) {
      data_at_GaussPt(gg,rr) = cblas_ddot(nb_dofs/rank,&data.getN()(gg,0),1,&data.getFieldData()[rr],rank);
      for(int dd = 0;dd<dim;dd++) {
	dataGrad_at_GaussPt(gg,dim*rr+dd) = cblas_ddot(nb_dofs/rank,&data.getDiffN()(gg,dd),1,&data.getFieldData()[rr],rank);
      }
    }
  }

  PetscFunctionReturn(0);
}

PetscErrorCode VolumeH1H1ElementForcesAndSurcesCore::operator()() {
  PetscFunctionBegin;

  ierr = getEdgesSense(data); CHKERRQ(ierr);
  ierr = getFacesSense(data); CHKERRQ(ierr);
  ierr = getEdgesOrder(data); CHKERRQ(ierr);
  ierr = getFacesOrder(data); CHKERRQ(ierr);
  ierr = getOrderVolume(data); CHKERRQ(ierr);
  ierr = getFaceNodes(data); CHKERRQ(ierr);

  int order = 1;
  for(int ee = 0;ee<data.eDges.size();ee++) {
    order = max(order,data.eDges[ee].getOrder());
  }

  int nb_gauss_pts = gm_rule_size(order,3);
  gaussPts.resize(4,nb_gauss_pts);
  ierr = Grundmann_Moeller_integration_points_3D_TET(
    order,&gaussPts(0,0),&gaussPts(1,0),&gaussPts(2,0),&gaussPts(3,0)); CHKERRQ(ierr);
  ierr = shapeTETFunctions_H1(data,
    &gaussPts(0,0),&gaussPts(1,0),&gaussPts(2,0),nb_gauss_pts); CHKERRQ(ierr);

  EntityHandle ent = fe_ptr->get_ent();
  int num_nodes;
  const EntityHandle* conn;
  rval = mField.get_moab().get_connectivity(ent,conn,num_nodes,true); CHKERR_PETSC(rval);
  coords.resize(num_nodes*3);
  rval = mField.get_moab().get_coords(conn,num_nodes,&*coords.data().begin()); CHKERR_PETSC(rval);
  vOlume = Shape_intVolumeMBTET(&*data.nOdes[0].getDiffN().data().begin(),&*coords.data().begin()); 
  invJac.resize(3,3);
  ierr = ShapeJacMBTET(&*data.nOdes[0].getDiffN().data().begin(),&*coords.begin(),&*invJac.data().begin()); CHKERRQ(ierr);
  ierr = Shape_invJac(&*invJac.data().begin()); CHKERRQ(ierr);

  coordsAtGaussPts.resize(nb_gauss_pts,3);
  for(int gg = 0;gg<nb_gauss_pts;gg++) {
    for(int dd = 0;dd<3;dd++) {
      coordsAtGaussPts(gg,dd) = cblas_ddot(4,&data.nOdes[0].getN()(gg,0),1,&coords[dd],3);
    }
  }

  try {
    ierr = opSetJac.opNH1(data); CHKERRQ(ierr);
  } catch (exception& ex) {
    ostringstream ss;
    ss << "thorw in method: " << ex.what() << " at line " << __LINE__ << " in file " << __FILE__;
    SETERRQ(PETSC_COMM_SELF,1,ss.str().c_str());
  }

  for(
    vector<UserDataOperator*>::iterator oit = vecUserOpNH1.begin();
    oit != vecUserOpNH1.end(); oit++) {

    (*oit)->ptrFE = this;

    ierr = getRowNodesIndices(data,(*oit)->row_field_name); CHKERRQ(ierr);
    ierr = getEdgeRowIndices(data,(*oit)->row_field_name); CHKERRQ(ierr);
    ierr = getFacesRowIndices(data,(*oit)->row_field_name); CHKERRQ(ierr);
    ierr = getTetRowIndices(data,(*oit)->row_field_name); CHKERRQ(ierr);

    ierr = getNodesFieldData(data,(*oit)->col_field_name); CHKERRQ(ierr);
    ierr = getEdgeFieldData(data,(*oit)->col_field_name); CHKERRQ(ierr);
    ierr = getFacesFieldData(data,(*oit)->col_field_name); CHKERRQ(ierr);
    ierr = getTetFieldData(data,(*oit)->col_field_name); CHKERRQ(ierr);

    try {
      ierr = (*oit)->opNH1(data); CHKERRQ(ierr);
    } catch (exception& ex) {
      ostringstream ss;
      ss << "thorw in method: " << ex.what() << " at line " << __LINE__ << " in file " << __FILE__;
      SETERRQ(PETSC_COMM_SELF,1,ss.str().c_str());
    }

  }

  for(
    vector<UserDataOperator*>::iterator oit = vecUserOpNH1NH1.begin();
    oit != vecUserOpNH1NH1.end(); oit++) {

    (*oit)->ptrFE = this;

    ierr = getRowNodesIndices(data,(*oit)->row_field_name); CHKERRQ(ierr);
    ierr = getEdgeRowIndices(data,(*oit)->row_field_name); CHKERRQ(ierr);
    ierr = getFacesRowIndices(data,(*oit)->row_field_name); CHKERRQ(ierr);
    ierr = getTetRowIndices(data,(*oit)->row_field_name); CHKERRQ(ierr);

    ierr = getColNodesIndices(derived_data,(*oit)->col_field_name); CHKERRQ(ierr);
    ierr = getEdgeColIndices(derived_data,(*oit)->col_field_name); CHKERRQ(ierr);
    ierr = getFacesColIndices(derived_data,(*oit)->col_field_name); CHKERRQ(ierr);
    ierr = getTetColIndices(derived_data,(*oit)->col_field_name); CHKERRQ(ierr);

    ierr = getNodesFieldData(data,(*oit)->col_field_name); CHKERRQ(ierr);
    ierr = getEdgeFieldData(data,(*oit)->col_field_name); CHKERRQ(ierr);
    ierr = getFacesFieldData(data,(*oit)->col_field_name); CHKERRQ(ierr);
    ierr = getTetFieldData(data,(*oit)->col_field_name); CHKERRQ(ierr);

    try {
      ierr = (*oit)->opNH1NH1(data,derived_data); CHKERRQ(ierr);
    } catch (exception& ex) {
      ostringstream ss;
      ss << "thorw in method: " << ex.what() << " at line " << __LINE__ << " in file " << __FILE__;
      SETERRQ(PETSC_COMM_SELF,1,ss.str().c_str());
    }

  }


  PetscFunctionReturn(0);
}

}
