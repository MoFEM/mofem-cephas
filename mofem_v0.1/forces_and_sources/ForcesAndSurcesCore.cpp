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
using namespace MoFEM;

PetscErrorCode ForcesAndSurcesCore::getEdgesSense(dataForcesAndSurcesCore &data) {
    PetscFunctionBegin;
    SideNumber_multiIndex& side_table = const_cast<SideNumber_multiIndex&>(fe_ptr->get_side_number_table());
    SideNumber_multiIndex::nth_index<2>::type::iterator siit = side_table.get<2>().lower_bound(MBEDGE);
    SideNumber_multiIndex::nth_index<2>::type::iterator hi_siit = side_table.get<2>().upper_bound(MBEDGE);
    data.edgesSense.resize(distance(siit,hi_siit));
    for(;siit!=hi_siit;siit++) {
      data.edgesSense[siit->side_number] = siit->sense;
    }
    PetscFunctionReturn(0);
  }

PetscErrorCode ForcesAndSurcesCore::getFacesSense(dataForcesAndSurcesCore &data) {
    PetscFunctionBegin;
    SideNumber_multiIndex& side_table = const_cast<SideNumber_multiIndex&>(fe_ptr->get_side_number_table());
    SideNumber_multiIndex::nth_index<2>::type::iterator siit = side_table.get<2>().lower_bound(MBTRI);
    SideNumber_multiIndex::nth_index<2>::type::iterator hi_siit = side_table.get<2>().upper_bound(MBTRI);
    data.facesSense.resize(distance(siit,hi_siit));
    for(;siit!=hi_siit;siit++) {
      data.facesSense[siit->side_number] = siit->sense;
    }
    PetscFunctionReturn(0);
  }

PetscErrorCode ForcesAndSurcesCore::getEdgesOrder(dataForcesAndSurcesCore &data,const string &field_name) {
    PetscFunctionBegin;
    SideNumber_multiIndex& side_table = const_cast<SideNumber_multiIndex&>(fe_ptr->get_side_number_table());
    data.edgesOrder.resize(side_table.get<2>().count(MBEDGE));
    FEDofMoFEMEntity_multiIndex::index<Composite_Name_And_Type_mi_tag>::type &data_dofs =
      const_cast<FEDofMoFEMEntity_multiIndex::index<Composite_Name_And_Type_mi_tag>::type&>(
	fe_ptr->get_data_dofs().get<Composite_Name_And_Type_mi_tag>());
    FEDofMoFEMEntity_multiIndex::index<Composite_Name_And_Type_mi_tag>::type::iterator dit,hi_dit;
    dit = data_dofs.lower_bound(boost::make_tuple(field_name,MBEDGE));
    hi_dit = data_dofs.upper_bound(boost::make_tuple(field_name,MBEDGE));
    for(;dit!=hi_dit;dit++) {
      ApproximationOrder ent_order = dit->get_max_order();
      int side_number = dit->side_number_ptr->side_number;
      data.edgesOrder[side_number] = data.edgesOrder[side_number] > ent_order ? data.edgesOrder[side_number] : ent_order;
    }
    PetscFunctionReturn(0);
  }

PetscErrorCode ForcesAndSurcesCore::getFacesOrder(dataForcesAndSurcesCore &data,const string &field_name) {
    PetscFunctionBegin;
    SideNumber_multiIndex& side_table = const_cast<SideNumber_multiIndex&>(fe_ptr->get_side_number_table());
    data.facesOrder.resize(side_table.get<2>().count(MBTRI));
    FEDofMoFEMEntity_multiIndex::index<Composite_Name_And_Type_mi_tag>::type &data_dofs =
      const_cast<FEDofMoFEMEntity_multiIndex::index<Composite_Name_And_Type_mi_tag>::type&>(
	fe_ptr->get_data_dofs().get<Composite_Name_And_Type_mi_tag>());
    FEDofMoFEMEntity_multiIndex::index<Composite_Name_And_Type_mi_tag>::type::iterator dit,hi_dit;
    dit = data_dofs.lower_bound(boost::make_tuple(field_name,MBTRI));
    hi_dit = data_dofs.upper_bound(boost::make_tuple(field_name,MBTRI));
    for(;dit!=hi_dit;dit++) {
      ApproximationOrder ent_order = dit->get_max_order();
      int side_number = dit->side_number_ptr->side_number;
      data.facesOrder[side_number] = data.facesOrder[side_number] > ent_order ? data.facesOrder[side_number] : ent_order;
    }
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
    data.volumeOrder = -1;
    for(;dit!=hi_dit;dit++) {
      ApproximationOrder ent_order = dit->get_max_order();
      int side_number = dit->side_number_ptr->side_number;
      data.volumeOrder = data.volumeOrder > ent_order ? data.volumeOrder : ent_order;
    }
    PetscFunctionReturn(0);
  }

PetscErrorCode ForcesAndSurcesCore::getNodesIndices(
    dataForcesAndSurcesCore &data,
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
      data,field_name,const_cast<FENumeredDofMoFEMEntity_multiIndex&>(fe_ptr->get_rows_dofs()),data.nodesIndices); CHKERRQ(ierr);
    PetscFunctionReturn(0);
  }

PetscErrorCode ForcesAndSurcesCore::getColNodesIndices(dataForcesAndSurcesCore &data,const string &field_name) {
    PetscFunctionBegin;
    PetscErrorCode ierr;
    ierr = getNodesIndices(data,field_name,const_cast<FENumeredDofMoFEMEntity_multiIndex&>(fe_ptr->get_cols_dofs()),data.nodesIndices); CHKERRQ(ierr);
    PetscFunctionReturn(0);
  }

PetscErrorCode ForcesAndSurcesCore::getTypeIndices(
  dataForcesAndSurcesCore &data,
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
  dataForcesAndSurcesCore &data,
  const string &field_name,FENumeredDofMoFEMEntity_multiIndex &dofs,EntityType type,
  ublas::vector<ublas::vector<int> > &indices) {
    PetscFunctionBegin;
    PetscErrorCode ierr;
    SideNumber_multiIndex& side_table = const_cast<SideNumber_multiIndex&>(fe_ptr->get_side_number_table());
    SideNumber_multiIndex::nth_index<2>::type::iterator siit = side_table.get<2>().lower_bound(type);
    SideNumber_multiIndex::nth_index<2>::type::iterator hi_siit = side_table.get<2>().upper_bound(type);
    indices.resize(side_table.get<2>().count(type));
    for(;siit!=hi_siit;siit++) {
      ierr = getTypeIndices(data,field_name,dofs,type,siit->side_number,indices[siit->side_number]); CHKERRQ(ierr);
    }
    PetscFunctionReturn(0);
  }

PetscErrorCode ForcesAndSurcesCore::getEdgeRowIndices(dataForcesAndSurcesCore &data,const string &field_name) {
    PetscFunctionBegin;
    PetscErrorCode ierr;
    ierr = getTypeIndices(
      data,field_name,const_cast<FENumeredDofMoFEMEntity_multiIndex&>(fe_ptr->get_rows_dofs()),MBEDGE,data.edgesIndcies); CHKERRQ(ierr);
    PetscFunctionReturn(0);
  }

PetscErrorCode ForcesAndSurcesCore::getEdgeColIndices(dataForcesAndSurcesCore &data,const string &field_name) {
    PetscFunctionBegin;
    PetscErrorCode ierr;
    ierr = getTypeIndices(
      data,field_name,const_cast<FENumeredDofMoFEMEntity_multiIndex&>(fe_ptr->get_cols_dofs()),MBEDGE,data.edgesIndcies); CHKERRQ(ierr);
    PetscFunctionReturn(0);
  }

PetscErrorCode ForcesAndSurcesCore::getFacesRowIndices(
    dataForcesAndSurcesCore &data,const string &field_name) {
    PetscFunctionBegin;
    PetscErrorCode ierr;
    ierr = getTypeIndices(data,field_name,const_cast<FENumeredDofMoFEMEntity_multiIndex&>(fe_ptr->get_rows_dofs()),MBTRI,data.facesIndices); CHKERRQ(ierr);
    PetscFunctionReturn(0);
  }

PetscErrorCode ForcesAndSurcesCore::getFacesColIndices(dataForcesAndSurcesCore &data,const string &field_name) {
    PetscFunctionBegin;
    PetscErrorCode ierr;
    ierr = getTypeIndices(data,field_name,const_cast<FENumeredDofMoFEMEntity_multiIndex&>(fe_ptr->get_cols_dofs()),MBTRI,data.facesIndices); CHKERRQ(ierr);
    PetscFunctionReturn(0);
  }

PetscErrorCode ForcesAndSurcesCore::getTetRowIndices(dataForcesAndSurcesCore &data,const string &field_name) {
    PetscFunctionBegin;
    PetscErrorCode ierr;
    ierr = getTypeIndices(data,field_name,const_cast<FENumeredDofMoFEMEntity_multiIndex&>(fe_ptr->get_rows_dofs()),MBTET,0,data.volumeIndices); CHKERRQ(ierr);
    PetscFunctionReturn(0);
  }

PetscErrorCode ForcesAndSurcesCore::getTetColIndices(dataForcesAndSurcesCore &data,const string &field_name) {
    PetscFunctionBegin;
    PetscErrorCode ierr;
    ierr = getTypeIndices(data,field_name,const_cast<FENumeredDofMoFEMEntity_multiIndex&>(fe_ptr->get_cols_dofs()),MBTET,0,data.volumeIndices); CHKERRQ(ierr);
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

    data.nodesNH1.resize(G_DIM,4);
    ierr = ShapeMBTET(&*data.nodesNH1.data().begin(),G_X,G_Y,G_Z,G_DIM); CHKERRQ(ierr);
    data.diffNodesNH1.resize(4,3);
    ierr = ShapeDiffMBTET(&*data.diffNodesNH1.data().begin()); CHKERRQ(ierr);

    //edges
    if(data.edgesSense.size()!=6) {
      SETERRQ(PETSC_COMM_SELF,1,"data inconsistency");
    }
    if(data.edgesOrder.size()!=6) {
      SETERRQ(PETSC_COMM_SELF,1,"data inconsistency");
    }
    data.edgesNH1.resize(6);
    double *_H1edgeN_[6];
    for(int ee = 0;ee<6;ee++) {
      data.edgesNH1[ee].resize(G_DIM,NBEDGE_H1(data.edgesOrder[ee]));
      _H1edgeN_[ee] = &*data.edgesNH1[ee].data().begin();
    }
    ierr = H1_EdgeShapeFunctions_MBTET(
      &*data.edgesSense.data().begin(),&*data.edgesOrder.data().begin(),&*data.nodesNH1.data().begin(),data.diffNodesNH1.data().begin(),_H1edgeN_,NULL,G_DIM); CHKERRQ(ierr);

    //faces
    if(data.facesOrder.size()!=4) {
      SETERRQ(PETSC_COMM_SELF,1,"data inconsistency");
    }
    data.facesNH1.resize(4);
    double *_H1faceN_[4];
    for(int ff = 0;ff<4;ff++) {
      data.facesNH1[ff].resize(G_DIM,NBFACE_H1(data.facesOrder[ff]));
      _H1faceN_[ff] = &*data.facesNH1[ff].data().begin();
    }
    if(data.facesNodes.size1() != 4) {
      SETERRQ(PETSC_COMM_SELF,1,"data inconsistency");
    }
    if(data.facesNodes.size2() != 3) {
      SETERRQ(PETSC_COMM_SELF,1,"data inconsistency");
    }
    ierr = H1_FaceShapeFunctions_MBTET(
      &*data.facesNodes.data().begin(),&*data.facesOrder.data().begin(),
      &*data.nodesNH1.data().begin(),data.diffNodesNH1.data().begin(),_H1faceN_,NULL,G_DIM); CHKERRQ(ierr);

    //volume
    data.volumeNH1.resize(G_DIM,NBVOLUME_H1(data.volumeOrder));
    ierr = H1_VolumeShapeFunctions_MBTET(
      data.volumeOrder,&*data.nodesNH1.data().begin(),data.diffNodesNH1.data().begin(),&*data.volumeNH1.data().begin(),NULL,G_DIM); CHKERRQ(ierr);

    PetscFunctionReturn(0);
  }
  
PetscErrorCode ForcesAndSurcesCore::shapeTETFunctions_L2(
    dataForcesAndSurcesCore &data,
    const double *G_X,const double *G_Y,const double *G_Z,const int G_DIM) {
    PetscFunctionBegin;

    PetscErrorCode ierr;

    data.nodesNL2.resize(G_DIM,4);
    ierr = ShapeMBTET(&*data.nodesNL2.data().begin(),G_X,G_Y,G_Z,G_DIM); CHKERRQ(ierr);
    data.diffNodesNL2.resize(4,3);
    ierr = ShapeDiffMBTET(&*data.diffNodesNL2.data().begin()); CHKERRQ(ierr);

    data.volumeNL2.resize(G_DIM,NBVOLUME_L2(data.volumeOrder));
    ierr = L2_ShapeFunctions_MBTET(
      data.volumeOrder,&*data.nodesNL2.data().begin(),&*data.diffNodesNL2.data().begin(),&*data.volumeNL2.data().begin(),NULL,G_DIM); CHKERRQ(ierr);
    PetscFunctionReturn(0);
  }

PetscErrorCode ForcesAndSurcesCore::shapeTRIFunctions_H1(
    dataForcesAndSurcesCore &data,
    ApproximationOrder order,const double *G_X,const double *G_Y,const int G_DIM) {
    PetscFunctionBegin;
    
    SETERRQ(PETSC_COMM_SELF,1,"not implemented");

    PetscFunctionReturn(0);
  }


PetscErrorCode dataOperator::operator()(
    dataForcesAndSurcesCore &row_data,dataForcesAndSurcesCore &col_data) {
  PetscFunctionBegin;
  PetscErrorCode ierr;

  int G_DIM = row_data.nodesNH1.size1();
  for(int gg = 0;gg<G_DIM;gg++) {
    ierr = doWork(
	gg,-1,-1,MBVERTEX,MBVERTEX,
	row_data.nodesIndices,col_data.nodesIndices,
	row_data.nodesNH1,col_data.nodesNH1); CHKERRQ(ierr);
  }

  for(int ee = 0;ee<row_data.edgesNH1.size();ee++) {
    int G_DIM = row_data.edgesNH1[ee].size1();
    int G_DIM_NODES = row_data.nodesNH1.size1();
    if(G_DIM != G_DIM_NODES) {
      SETERRQ(PETSC_COMM_SELF,1,"data inconsistency");
    }
    for(int gg = 0;gg<G_DIM;gg++) {
      ierr = doWork(
	gg,ee,-1,MBEDGE,MBVERTEX,
	row_data.edgesIndcies[ee],col_data.nodesIndices,
	row_data.edgesNH1[ee],col_data.nodesNH1); CHKERRQ(ierr);
    }
    int G_DIM_VOLUME = row_data.volumeNH1.size1();
    if(G_DIM != G_DIM_VOLUME) {
      SETERRQ(PETSC_COMM_SELF,1,"data inconsistency");
    }
    for(int gg = 0;gg<G_DIM;gg++) {
      ierr = doWork(
	gg,ee,-1,MBEDGE,MBTET,
	row_data.edgesIndcies[ee],col_data.volumeIndices,
	row_data.edgesNH1[ee],col_data.volumeNH1); CHKERRQ(ierr);
    }
    for(int EE = 0;EE<row_data.edgesNH1.size();EE++) {
      int G_DIM_EE = row_data.edgesNH1[EE].size1();
      if(G_DIM != G_DIM_EE) {
	SETERRQ(PETSC_COMM_SELF,1,"data inconsistency");
      }
      for(int gg = 0;gg<G_DIM;gg++) {
	ierr = doWork(
	  gg,ee,EE,MBEDGE,MBEDGE,
	  row_data.edgesIndcies[ee],col_data.edgesIndcies[EE],
	  row_data.edgesNH1[ee],col_data.edgesNH1[EE]); CHKERRQ(ierr);
      }
    }
    for(int FF = 0;FF<row_data.facesNH1.size();FF++) {
      int G_DIM_FF = row_data.facesNH1[FF].size1();
      if(G_DIM != G_DIM_FF) {
	SETERRQ(PETSC_COMM_SELF,1,"data inconsistency");
      }
      for(int gg = 0;gg<G_DIM;gg++) {
	ierr = doWork(
	  gg,ee,FF,MBEDGE,MBTRI,
	  row_data.edgesIndcies[ee],col_data.facesIndices[FF],
	  row_data.edgesNH1[ee],col_data.facesNH1[FF]); CHKERRQ(ierr);
      }
    }
  }

  for(int ff = 0;ff<row_data.facesNH1.size();ff++) {
    int G_DIM = row_data.facesNH1[ff].size1();
    int G_DIM_NODES = row_data.nodesNH1.size1();
    if(G_DIM != G_DIM_NODES) {
      SETERRQ(PETSC_COMM_SELF,1,"data inconsistency");
    }
    for(int gg = 0;gg<G_DIM;gg++) {
      ierr = doWork(
	gg,ff,-1,MBTRI,MBVERTEX,
	row_data.facesIndices[ff],col_data.nodesIndices,
	row_data.facesNH1[ff],col_data.nodesNH1); CHKERRQ(ierr);
    }
    int G_DIM_VOLUME = row_data.volumeNH1.size1();
    if(G_DIM != G_DIM_VOLUME) {
      SETERRQ(PETSC_COMM_SELF,1,"data inconsistency");
    }
    for(int gg = 0;gg<G_DIM;gg++) {
      ierr = doWork(
	gg,ff,-1,MBTRI,MBTET,
	row_data.facesIndices[ff],col_data.volumeIndices,
	row_data.facesNH1[ff],col_data.volumeNH1); CHKERRQ(ierr);
    }
    for(int EE = 0;EE<row_data.edgesNH1.size();EE++) {
      int G_DIM_EE = row_data.edgesNH1[EE].size1();
      if(G_DIM != G_DIM_EE) {
	SETERRQ(PETSC_COMM_SELF,1,"data inconsistency");
      }
      for(int gg = 0;gg<G_DIM;gg++) {
	ierr = doWork(
	  gg,ff,EE,MBTRI,MBEDGE,
	  row_data.facesIndices[ff],col_data.edgesIndcies[EE],
	  row_data.facesNH1[ff],col_data.edgesNH1[EE]); CHKERRQ(ierr);
      }
    }
    for(int FF = 0;FF<row_data.facesNH1.size();FF++) {
      int G_DIM_FF = row_data.facesNH1[ff].size1();
      if(G_DIM != G_DIM_FF) {
	SETERRQ(PETSC_COMM_SELF,1,"data inconsistency");
      }
      for(int gg = 0;gg<G_DIM;gg++) {
	ierr = doWork(
	  gg,ff,FF,MBTRI,MBTRI,
	  row_data.facesIndices[ff],col_data.facesIndices[FF],
	  row_data.facesNH1[ff],col_data.facesNH1[FF]); CHKERRQ(ierr);
      }
    }
  }


  {

    int G_DIM = row_data.volumeNH1.size1();
    for(int gg = 0;gg<G_DIM;gg++) {
      ierr = doWork(
	gg,-1,-1,MBTET,MBTET,
	row_data.volumeIndices,col_data.volumeIndices,
	row_data.volumeNH1,col_data.volumeNH1); CHKERRQ(ierr);
    }
    int nb_cols_nodes = row_data.nodesNH1.size2();
    for(int gg = 0;gg<G_DIM;gg++) {
      ierr = doWork(
	gg,-1,-1,MBTET,MBVERTEX,
	row_data.volumeIndices,col_data.nodesIndices,
	row_data.volumeNH1,col_data.nodesNH1); CHKERRQ(ierr);
    }

  }


  PetscFunctionReturn(0);
}

PetscErrorCode mult_H1_H1::doWork(
    int gg,int side1,int side2,
    EntityType row_type,EntityType col_type,
    ublas::vector<DofIdx> &row_indices,ublas::vector<DofIdx> &col_indices,
    ublas::matrix<FieldData> &rows_N,ublas::matrix<FieldData> &cols_N) {
  PetscFunctionBegin;

  int nb_row_dofs = rows_N.size2();
  int nb_col_dofs = rows_N.size2();

  NN.resize(nb_row_dofs,nb_col_dofs);
  bzero(NN.data().begin(),nb_row_dofs*nb_col_dofs*sizeof(FieldData));

  cblas_dger(CblasRowMajor,
    nb_row_dofs,nb_col_dofs,
      1,&rows_N(gg,0),1,&cols_N(gg,0),1,
      &*NN.data().begin(),nb_row_dofs);

  PetscFunctionReturn(0);


}





