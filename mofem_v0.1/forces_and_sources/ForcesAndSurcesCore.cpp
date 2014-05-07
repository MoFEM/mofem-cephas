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

using namespace boost::numeric;

namespace MoFEM {

ostream& operator<<(ostream& os,const dataForcesAndSurcesCore::entData &e) {
  os << 
    "sEnse: " << e.sEnse << endl << 
    "oRder: " << e.oRder << endl <<
    "iNdices: " << e.iNdices << endl;
  os.precision(2);
  os << 
    "fieldData: " << e.fieldData << endl;
  os <<
    "N: " << e.N << endl <<
    "diffN: " << e.diffN;
  return os;
}

ostream& operator<<(ostream& os,const dataForcesAndSurcesCore &e) {
  os << 
    "nOdes: " << endl << e.nOdes << endl;
  for(int ee = 0;ee < e.eDges.size(); ee++) {
    os << "eDges[" << ee << "]" << endl << e.eDges[ee] << endl;
  }
  for(int ff = 0;ff < e.fAces.size(); ff++) {
    os << "fAces[" << ff << "] " << endl << e.fAces[ff] << endl;
  }
  os << 
    "vOlume: " << endl << e.vOlume;
  return os;
}

PetscErrorCode ForcesAndSurcesCore::getSense(EntityType type,vector<dataForcesAndSurcesCore::entData> &data) {
    PetscFunctionBegin;
    SideNumber_multiIndex& side_table = const_cast<SideNumber_multiIndex&>(fe_ptr->get_side_number_table());
    SideNumber_multiIndex::nth_index<2>::type::iterator siit = side_table.get<2>().lower_bound(type);
    SideNumber_multiIndex::nth_index<2>::type::iterator hi_siit = side_table.get<2>().upper_bound(type);
    data.resize(distance(siit,hi_siit));
    for(;siit!=hi_siit;siit++) {
      data[siit->side_number].sEnse = siit->sense;
    }
    PetscFunctionReturn(0);
  }

PetscErrorCode ForcesAndSurcesCore::getEdgesSense(dataForcesAndSurcesCore &data) {
    PetscFunctionBegin;
    PetscErrorCode ierr;
    //PetscAttachDebugger();
    ierr = getSense(MBEDGE,data.eDges); CHKERRQ(ierr);
    PetscFunctionReturn(0);
  }

PetscErrorCode ForcesAndSurcesCore::getFacesSense(dataForcesAndSurcesCore &data) {
    PetscFunctionBegin;
    PetscErrorCode ierr;
    ierr = getSense(MBTRI,data.fAces); CHKERRQ(ierr);
    PetscFunctionReturn(0);
  }

PetscErrorCode ForcesAndSurcesCore::getOrder(EntityType type,vector<dataForcesAndSurcesCore::entData> &data,const string &field_name) {
    PetscFunctionBegin;
    SideNumber_multiIndex& side_table = const_cast<SideNumber_multiIndex&>(fe_ptr->get_side_number_table());
    data.resize(side_table.get<2>().count(type));
    FEDofMoFEMEntity_multiIndex::index<Composite_Name_And_Type_mi_tag>::type &data_dofs =
      const_cast<FEDofMoFEMEntity_multiIndex::index<Composite_Name_And_Type_mi_tag>::type&>(
	fe_ptr->get_data_dofs().get<Composite_Name_And_Type_mi_tag>());
    FEDofMoFEMEntity_multiIndex::index<Composite_Name_And_Type_mi_tag>::type::iterator dit,hi_dit;
    dit = data_dofs.lower_bound(boost::make_tuple(field_name,type));
    hi_dit = data_dofs.upper_bound(boost::make_tuple(field_name,type));
    for(;dit!=hi_dit;dit++) {
      ApproximationOrder ent_order = dit->get_max_order();
      int side_number = dit->side_number_ptr->side_number;
      data[side_number].oRder = data[side_number].oRder > ent_order ? data[side_number].oRder : ent_order;
    }
    PetscFunctionReturn(0);
  }


PetscErrorCode ForcesAndSurcesCore::getEdgesOrder(dataForcesAndSurcesCore &data,const string &field_name) {
    PetscFunctionBegin;
    PetscErrorCode ierr;
    ierr = getOrder(MBEDGE,data.eDges,field_name); CHKERRQ(ierr);
    PetscFunctionReturn(0);
  }

PetscErrorCode ForcesAndSurcesCore::getFacesOrder(dataForcesAndSurcesCore &data,const string &field_name) {
    PetscFunctionBegin;
    PetscErrorCode ierr;
    ierr = getOrder(MBTRI,data.fAces,field_name); CHKERRQ(ierr);
    PetscFunctionReturn(0);
  }

PetscErrorCode ForcesAndSurcesCore::getOrderVolume(dataForcesAndSurcesCore &data,const string &field_name) {
    PetscFunctionBegin;
    FEDofMoFEMEntity_multiIndex::index<Composite_Name_And_Type_mi_tag>::type &data_dofs =
      const_cast<FEDofMoFEMEntity_multiIndex::index<Composite_Name_And_Type_mi_tag>::type&>(
	fe_ptr->get_data_dofs().get<Composite_Name_And_Type_mi_tag>());
    FEDofMoFEMEntity_multiIndex::index<Composite_Name_And_Type_mi_tag>::type::iterator dit,hi_dit;
    dit = data_dofs.lower_bound(boost::make_tuple(field_name,MBTET));
    hi_dit = data_dofs.upper_bound(boost::make_tuple(field_name,MBTET));
    data.vOlume.oRder = -1;
    for(;dit!=hi_dit;dit++) {
      ApproximationOrder ent_order = dit->get_max_order();
      int side_number = dit->side_number_ptr->side_number;
      data.vOlume.oRder = data.vOlume.oRder > ent_order ? data.vOlume.oRder : ent_order;
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

PetscErrorCode ForcesAndSurcesCore::getRowNodesIndices(dataForcesAndSurcesCore &data,const string &field_name) {
    PetscFunctionBegin;
    PetscErrorCode ierr;
    ierr = getNodesIndices(
      field_name,const_cast<FENumeredDofMoFEMEntity_multiIndex&>(fe_ptr->get_rows_dofs()),data.nOdes.iNdices); CHKERRQ(ierr);
    PetscFunctionReturn(0);
  }

PetscErrorCode ForcesAndSurcesCore::getColNodesIndices(dataForcesAndSurcesCore &data,const string &field_name) {
    PetscFunctionBegin;
    PetscErrorCode ierr;
    ierr = getNodesIndices(field_name,const_cast<FENumeredDofMoFEMEntity_multiIndex&>(fe_ptr->get_cols_dofs()),data.nOdes.iNdices); CHKERRQ(ierr);
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
  vector<dataForcesAndSurcesCore::entData> &data) {
    PetscFunctionBegin;
    PetscErrorCode ierr;
    SideNumber_multiIndex& side_table = const_cast<SideNumber_multiIndex&>(fe_ptr->get_side_number_table());
    SideNumber_multiIndex::nth_index<2>::type::iterator siit = side_table.get<2>().lower_bound(type);
    SideNumber_multiIndex::nth_index<2>::type::iterator hi_siit = side_table.get<2>().upper_bound(type);
    data.resize(side_table.get<2>().count(type));
    for(;siit!=hi_siit;siit++) {
      ierr = getTypeIndices(field_name,dofs,type,siit->side_number,data[siit->side_number].iNdices); CHKERRQ(ierr);
    }
    PetscFunctionReturn(0);
  }

PetscErrorCode ForcesAndSurcesCore::getEdgeRowIndices(dataForcesAndSurcesCore &data,const string &field_name) {
    PetscFunctionBegin;
    PetscErrorCode ierr;
    ierr = getTypeIndices(
      field_name,const_cast<FENumeredDofMoFEMEntity_multiIndex&>(fe_ptr->get_rows_dofs()),MBEDGE,data.eDges); CHKERRQ(ierr);
    PetscFunctionReturn(0);
  }

PetscErrorCode ForcesAndSurcesCore::getEdgeColIndices(dataForcesAndSurcesCore &data,const string &field_name) {
    PetscFunctionBegin;
    PetscErrorCode ierr;
    ierr = getTypeIndices(
      field_name,const_cast<FENumeredDofMoFEMEntity_multiIndex&>(fe_ptr->get_cols_dofs()),MBEDGE,data.eDges); CHKERRQ(ierr);
    PetscFunctionReturn(0);
  }

PetscErrorCode ForcesAndSurcesCore::getFacesRowIndices(
    dataForcesAndSurcesCore &data,const string &field_name) {
    PetscFunctionBegin;
    PetscErrorCode ierr;
    ierr = getTypeIndices(field_name,const_cast<FENumeredDofMoFEMEntity_multiIndex&>(fe_ptr->get_rows_dofs()),MBTRI,data.fAces); CHKERRQ(ierr);
    PetscFunctionReturn(0);
  }

PetscErrorCode ForcesAndSurcesCore::getFacesColIndices(dataForcesAndSurcesCore &data,const string &field_name) {
    PetscFunctionBegin;
    PetscErrorCode ierr;
    ierr = getTypeIndices(field_name,const_cast<FENumeredDofMoFEMEntity_multiIndex&>(fe_ptr->get_cols_dofs()),MBTRI,data.fAces); CHKERRQ(ierr);
    PetscFunctionReturn(0);
  }

PetscErrorCode ForcesAndSurcesCore::getTetRowIndices(dataForcesAndSurcesCore &data,const string &field_name) {
    PetscFunctionBegin;
    PetscErrorCode ierr;
    ierr = getTypeIndices(field_name,const_cast<FENumeredDofMoFEMEntity_multiIndex&>(fe_ptr->get_rows_dofs()),MBTET,0,data.vOlume.iNdices); CHKERRQ(ierr);
    PetscFunctionReturn(0);
  }

PetscErrorCode ForcesAndSurcesCore::getTetColIndices(dataForcesAndSurcesCore &data,const string &field_name) {
    PetscFunctionBegin;
    PetscErrorCode ierr;
    ierr = getTypeIndices(field_name,const_cast<FENumeredDofMoFEMEntity_multiIndex&>(fe_ptr->get_cols_dofs()),MBTET,0,data.vOlume.iNdices); CHKERRQ(ierr);
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

PetscErrorCode ForcesAndSurcesCore::getNodesFieldData(dataForcesAndSurcesCore &data,const string &field_name) {
    PetscFunctionBegin;
    PetscErrorCode ierr;
    ierr = getNodesFieldData(
      field_name,const_cast<FEDofMoFEMEntity_multiIndex&>(fe_ptr->get_data_dofs()),data.nOdes.fieldData); CHKERRQ(ierr);
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
  vector<dataForcesAndSurcesCore::entData> &data) {
    PetscFunctionBegin;
    PetscErrorCode ierr;
    SideNumber_multiIndex& side_table = const_cast<SideNumber_multiIndex&>(fe_ptr->get_side_number_table());
    SideNumber_multiIndex::nth_index<2>::type::iterator siit = side_table.get<2>().lower_bound(type);
    SideNumber_multiIndex::nth_index<2>::type::iterator hi_siit = side_table.get<2>().upper_bound(type);
    data.resize(side_table.get<2>().count(type));
    for(;siit!=hi_siit;siit++) {
      ierr = getTypeFieldData(field_name,dofs,type,siit->side_number,data[siit->side_number].fieldData); CHKERRQ(ierr);
    }
    PetscFunctionReturn(0);
  }

PetscErrorCode ForcesAndSurcesCore::getEdgeFieldData(dataForcesAndSurcesCore &data,const string &field_name) {
    PetscFunctionBegin;
    PetscErrorCode ierr;
    ierr = getTypeFieldData(
      field_name,const_cast<FEDofMoFEMEntity_multiIndex&>(fe_ptr->get_data_dofs()),MBEDGE,data.eDges); CHKERRQ(ierr);
    PetscFunctionReturn(0);
  }

PetscErrorCode ForcesAndSurcesCore::getFacesFieldData(
    dataForcesAndSurcesCore &data,const string &field_name) {
    PetscFunctionBegin;
    PetscErrorCode ierr;
    ierr = getTypeFieldData(field_name,const_cast<FEDofMoFEMEntity_multiIndex&>(fe_ptr->get_data_dofs()),MBTRI,data.fAces); CHKERRQ(ierr);
    PetscFunctionReturn(0);
  }

PetscErrorCode ForcesAndSurcesCore::getTetFieldData(dataForcesAndSurcesCore &data,const string &field_name) {
    PetscFunctionBegin;
    PetscErrorCode ierr;
    ierr = getTypeFieldData(field_name,const_cast<FEDofMoFEMEntity_multiIndex&>(fe_ptr->get_data_dofs()),MBTET,0,data.vOlume.fieldData); CHKERRQ(ierr);
    PetscFunctionReturn(0);
  }

PetscErrorCode ForcesAndSurcesCore::getFaceNodes(dataForcesAndSurcesCore &data) {
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
    dataForcesAndSurcesCore &data,
    const double *G_X,const double *G_Y,const double *G_Z,const int G_DIM) {
    PetscFunctionBegin;

    PetscErrorCode ierr;

    data.nOdes.N.resize(G_DIM,4);
    ierr = ShapeMBTET(&*data.nOdes.N.data().begin(),G_X,G_Y,G_Z,G_DIM); CHKERRQ(ierr);
    data.nOdes.diffN.resize(4,3);
    ierr = ShapeDiffMBTET(&*data.nOdes.diffN.data().begin()); CHKERRQ(ierr);

    //edges
    if(data.eDges.size()!=6) {
      SETERRQ(PETSC_COMM_SELF,1,"data inconsistency");
    }
    if(data.eDges.size()!=6) {
      SETERRQ(PETSC_COMM_SELF,1,"data inconsistency");
    }
    int _sense_[6],_order_[6];
    double *_H1edgeN_[6],*_diffH1edgeN_[6];
    for(int ee = 0;ee<6;ee++) {
      if(data.eDges[ee].sEnse == 0) {
	SETERRQ(PETSC_COMM_SELF,1,"data inconsistency");
      }
      _sense_[ee] = data.eDges[ee].sEnse;
      _order_[ee] = data.eDges[ee].oRder;
      int nb_dofs = NBEDGE_H1(data.eDges[ee].oRder);
      data.eDges[ee].N.resize(G_DIM,nb_dofs);
      data.eDges[ee].diffN.resize(G_DIM,3*nb_dofs);
      _H1edgeN_[ee] = &*data.eDges[ee].N.data().begin();
      _diffH1edgeN_[ee] = &*data.eDges[ee].diffN.data().begin();
    }
    ierr = H1_EdgeShapeFunctions_MBTET(
      _sense_,_order_,
      &*data.nOdes.N.data().begin(),data.nOdes.diffN.data().begin(),
      _H1edgeN_,_diffH1edgeN_,G_DIM); CHKERRQ(ierr);

    //faces
    if(data.fAces.size()!=4) {
      SETERRQ(PETSC_COMM_SELF,1,"data inconsistency");
    }
    double *_H1faceN_[4],*_diffH1faceN_[4];
    for(int ff = 0;ff<4;ff++) {
      if(data.fAces[ff].sEnse == 0) {
	SETERRQ(PETSC_COMM_SELF,1,"data inconsistency");
      }
      int nb_dofs = NBFACE_H1(data.fAces[ff].oRder);
      _order_[ff] = data.fAces[ff].oRder;
      data.fAces[ff].N.resize(G_DIM,nb_dofs);
      data.fAces[ff].diffN.resize(G_DIM,3*nb_dofs);
      _H1faceN_[ff] = &*data.fAces[ff].N.data().begin();
      _diffH1faceN_[ff] = &*data.fAces[ff].diffN.data().begin();
    }
    if(data.facesNodes.size1() != 4) {
      SETERRQ(PETSC_COMM_SELF,1,"data inconsistency");
    }
    if(data.facesNodes.size2() != 3) {
      SETERRQ(PETSC_COMM_SELF,1,"data inconsistency");
    }
    ierr = H1_FaceShapeFunctions_MBTET(
      &*data.facesNodes.data().begin(),_order_,
      &*data.nOdes.N.data().begin(),data.nOdes.diffN.data().begin(),
      _H1faceN_,_diffH1faceN_,G_DIM); CHKERRQ(ierr);

    //volume
    int nb_vol_dofs = NBVOLUME_H1(data.vOlume.oRder);
    data.vOlume.N.resize(G_DIM,nb_vol_dofs);
    data.vOlume.diffN.resize(G_DIM,3*nb_vol_dofs);
    ierr = H1_VolumeShapeFunctions_MBTET(
      data.vOlume.oRder,&*data.nOdes.N.data().begin(),data.nOdes.diffN.data().begin(),
      &*data.vOlume.N.data().begin(),&*data.vOlume.diffN.data().begin(),G_DIM); CHKERRQ(ierr);

    PetscFunctionReturn(0);
  }
  
PetscErrorCode ForcesAndSurcesCore::shapeTETFunctions_L2(
    dataForcesAndSurcesCore &data,
    const double *G_X,const double *G_Y,const double *G_Z,const int G_DIM) {
    PetscFunctionBegin;

    PetscErrorCode ierr;

    data.nOdes.N.resize(G_DIM,4);
    ierr = ShapeMBTET(&*data.nOdes.N.data().begin(),G_X,G_Y,G_Z,G_DIM); CHKERRQ(ierr);
    data.nOdes.diffN.resize(4,3);
    ierr = ShapeDiffMBTET(&*data.nOdes.diffN.data().begin()); CHKERRQ(ierr);

    data.vOlume.N.resize(G_DIM,NBVOLUME_L2(data.vOlume.oRder));
    ierr = L2_ShapeFunctions_MBTET(
      data.vOlume.oRder,
      &*data.nOdes.N.data().begin(),&*data.nOdes.diffN.data().begin(),
      &*data.vOlume.N.data().begin(),NULL,G_DIM); CHKERRQ(ierr);
    PetscFunctionReturn(0);
  }

PetscErrorCode ForcesAndSurcesCore::shapeTRIFunctions_H1(
    dataForcesAndSurcesCore &data,
    ApproximationOrder order,const double *G_X,const double *G_Y,const int G_DIM) {
    PetscFunctionBegin;
    
    SETERRQ(PETSC_COMM_SELF,1,"not implemented");

    PetscFunctionReturn(0);
  }


PetscErrorCode dataOperator::opNH1NH1(
    dataForcesAndSurcesCore &row_data,dataForcesAndSurcesCore &col_data) {
  PetscFunctionBegin;
  PetscErrorCode ierr;

  int G_DIM = row_data.nOdes.N.size1();
  ierr = doWork(
    -1,-1,MBVERTEX,MBVERTEX,
    row_data.nOdes,col_data.nOdes); CHKERRQ(ierr);

  for(int ee = 0;ee<row_data.eDges.size();ee++) {
    int G_DIM = row_data.eDges[ee].N.size1();
    int G_DIM_NODES = col_data.nOdes.N.size1();
    if(G_DIM != G_DIM_NODES) {
      SETERRQ(PETSC_COMM_SELF,1,"data inconsistency");
    }
    ierr = doWork(
	ee,-1,MBEDGE,MBVERTEX,
	row_data.eDges[ee],col_data.nOdes); CHKERRQ(ierr);
    int G_DIM_VOLUME = col_data.vOlume.N.size1();
    if(G_DIM != G_DIM_VOLUME) {
      SETERRQ(PETSC_COMM_SELF,1,"data inconsistency");
    }
    ierr = doWork(
	ee,-1,MBEDGE,MBTET,
	row_data.eDges[ee],col_data.vOlume); CHKERRQ(ierr);
    for(int EE = 0;EE<col_data.eDges.size();EE++) {
      int G_DIM_EE = col_data.eDges[EE].N.size1();
      if(G_DIM != G_DIM_EE) {
	SETERRQ(PETSC_COMM_SELF,1,"data inconsistency");
      }
      ierr = doWork(
	  ee,EE,MBEDGE,MBEDGE,
	  row_data.eDges[ee],col_data.eDges[EE]); CHKERRQ(ierr);
    }
    for(int FF = 0;FF<col_data.fAces.size();FF++) {
      int G_DIM_FF = col_data.fAces[FF].N.size1();
      if(G_DIM != G_DIM_FF) {
	SETERRQ(PETSC_COMM_SELF,1,"data inconsistency");
      }
      ierr = doWork(
	  ee,FF,MBEDGE,MBTRI,
	  row_data.eDges[ee],col_data.fAces[FF]); CHKERRQ(ierr);
    }
  }

  for(int ff = 0;ff<row_data.fAces.size();ff++) {
    int G_DIM = row_data.fAces[ff].N.size1();
    int G_DIM_NODES = col_data.nOdes.N.size1();
    if(G_DIM != G_DIM_NODES) {
      SETERRQ(PETSC_COMM_SELF,1,"data inconsistency");
    }
    ierr = doWork(
	ff,-1,MBTRI,MBVERTEX,
	row_data.fAces[ff],col_data.nOdes); CHKERRQ(ierr);
    int G_DIM_VOLUME = col_data.vOlume.N.size1();
    if(G_DIM != G_DIM_VOLUME) {
      SETERRQ(PETSC_COMM_SELF,1,"data inconsistency");
    }
    ierr = doWork(
	ff,-1,MBTRI,MBTET,
	row_data.fAces[ff],col_data.vOlume); CHKERRQ(ierr);
    for(int EE = 0;EE<col_data.eDges.size();EE++) {
      int G_DIM_EE = col_data.eDges[EE].N.size1();
      if(G_DIM != G_DIM_EE) {
	SETERRQ(PETSC_COMM_SELF,1,"data inconsistency");
      }
      ierr = doWork(
	  ff,EE,MBTRI,MBEDGE,
	  row_data.fAces[ff],col_data.eDges[EE]); CHKERRQ(ierr);
    }
    for(int FF = 0;FF<col_data.fAces.size();FF++) {
      int G_DIM_FF = col_data.fAces[FF].N.size1();
      if(G_DIM != G_DIM_FF) {
	SETERRQ(PETSC_COMM_SELF,1,"data inconsistency");
      }
      ierr = doWork(
	  ff,FF,MBTRI,MBTRI,
	  row_data.fAces[ff],col_data.fAces[FF]); CHKERRQ(ierr);
    }
  }


  {

    int G_DIM = row_data.vOlume.N.size1();
    ierr = doWork(
	-1,-1,MBTET,MBTET,
	row_data.vOlume,col_data.vOlume); CHKERRQ(ierr);

  }

  PetscFunctionReturn(0);
}

PetscErrorCode dataOperator::opNH1(dataForcesAndSurcesCore &data) {
  PetscFunctionBegin;
  PetscErrorCode ierr;

  ierr = doWork(-1,MBVERTEX,data.nOdes); CHKERRQ(ierr);
  for(int ee = 0;ee<data.eDges.size();ee++) {
    ierr = doWork(ee,MBEDGE,data.eDges[ee]); CHKERRQ(ierr);
  }
  for(int ff = 0;ff<data.fAces.size();ff++) {
    ierr = doWork(ff,MBTRI,data.fAces[ff]); CHKERRQ(ierr);
  }
  ierr = doWork(-1,MBTET,data.vOlume); CHKERRQ(ierr);

  PetscFunctionReturn(0);
}

}
