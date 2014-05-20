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

PetscErrorCode ForcesAndSurcesCore::getEdgesSense() {
    PetscFunctionBegin;
    SideNumber_multiIndex& side_table = const_cast<SideNumber_multiIndex&>(fe_ptr->get_side_number_table());
    SideNumber_multiIndex::nth_index<2>::type::iterator siit = side_table.get<2>().lower_bound(MBEDGE);
    SideNumber_multiIndex::nth_index<2>::type::iterator hi_siit = side_table.get<2>().upper_bound(MBEDGE);
    edgesSense.resize(distance(siit,hi_siit));
    for(;siit!=hi_siit;siit++) {
      edgesSense[siit->side_number] = siit->sense;
    }
    PetscFunctionReturn(0);
  }

PetscErrorCode ForcesAndSurcesCore::getFacesSense() {
    PetscFunctionBegin;
    SideNumber_multiIndex& side_table = const_cast<SideNumber_multiIndex&>(fe_ptr->get_side_number_table());
    SideNumber_multiIndex::nth_index<2>::type::iterator siit = side_table.get<2>().lower_bound(MBTRI);
    SideNumber_multiIndex::nth_index<2>::type::iterator hi_siit = side_table.get<2>().upper_bound(MBTRI);
    facesSense.resize(distance(siit,hi_siit));
    for(;siit!=hi_siit;siit++) {
      facesSense[siit->side_number] = siit->sense;
    }
    PetscFunctionReturn(0);
  }

PetscErrorCode ForcesAndSurcesCore::getEdgesOrder(const string &field_name) {
    PetscFunctionBegin;
    SideNumber_multiIndex& side_table = const_cast<SideNumber_multiIndex&>(fe_ptr->get_side_number_table());
    edgesOrder.resize(side_table.get<2>().count(MBEDGE));
    FEDofMoFEMEntity_multiIndex::index<Composite_Name_And_Type_mi_tag>::type &data_dofs =
      const_cast<FEDofMoFEMEntity_multiIndex::index<Composite_Name_And_Type_mi_tag>::type&>(
	fe_ptr->get_data_dofs().get<Composite_Name_And_Type_mi_tag>());
    FEDofMoFEMEntity_multiIndex::index<Composite_Name_And_Type_mi_tag>::type::iterator dit,hi_dit;
    dit = data_dofs.lower_bound(boost::make_tuple(field_name,MBEDGE));
    hi_dit = data_dofs.upper_bound(boost::make_tuple(field_name,MBEDGE));
    for(;dit!=hi_dit;dit++) {
      ApproximationOrder ent_order = dit->get_max_order();
      int side_number = dit->side_number_ptr->side_number;
      edgesOrder[side_number] = edgesOrder[side_number] > ent_order ? edgesOrder[side_number] : ent_order;
    }
    PetscFunctionReturn(0);
  }

PetscErrorCode ForcesAndSurcesCore::getFacesOrder(const string &field_name) {
    PetscFunctionBegin;
    SideNumber_multiIndex& side_table = const_cast<SideNumber_multiIndex&>(fe_ptr->get_side_number_table());
    facesOrder.resize(side_table.get<2>().count(MBTRI));
    FEDofMoFEMEntity_multiIndex::index<Composite_Name_And_Type_mi_tag>::type &data_dofs =
      const_cast<FEDofMoFEMEntity_multiIndex::index<Composite_Name_And_Type_mi_tag>::type&>(
	fe_ptr->get_data_dofs().get<Composite_Name_And_Type_mi_tag>());
    FEDofMoFEMEntity_multiIndex::index<Composite_Name_And_Type_mi_tag>::type::iterator dit,hi_dit;
    dit = data_dofs.lower_bound(boost::make_tuple(field_name,MBTRI));
    hi_dit = data_dofs.upper_bound(boost::make_tuple(field_name,MBTRI));
    for(;dit!=hi_dit;dit++) {
      ApproximationOrder ent_order = dit->get_max_order();
      int side_number = dit->side_number_ptr->side_number;
      facesOrder[side_number] = facesOrder[side_number] > ent_order ? facesOrder[side_number] : ent_order;
    }
    PetscFunctionReturn(0);
  }

PetscErrorCode ForcesAndSurcesCore::getOrderVolume(const string &field_name) {
    PetscFunctionBegin;
    FEDofMoFEMEntity_multiIndex::index<Composite_Name_And_Type_mi_tag>::type &data_dofs =
      const_cast<FEDofMoFEMEntity_multiIndex::index<Composite_Name_And_Type_mi_tag>::type&>(
	fe_ptr->get_data_dofs().get<Composite_Name_And_Type_mi_tag>());
    FEDofMoFEMEntity_multiIndex::index<Composite_Name_And_Type_mi_tag>::type::iterator dit,hi_dit;
    dit = data_dofs.lower_bound(boost::make_tuple(field_name,MBTET));
    hi_dit = data_dofs.upper_bound(boost::make_tuple(field_name,MBTET));
    volumeOrder = -1;
    for(;dit!=hi_dit;dit++) {
      ApproximationOrder ent_order = dit->get_max_order();
      int side_number = dit->side_number_ptr->side_number;
      volumeOrder = volumeOrder > ent_order ? volumeOrder : ent_order;
    }
    PetscFunctionReturn(0);
  }

PetscErrorCode ForcesAndSurcesCore::getNodesIndices(const string &field_name,FENumeredDofMoFEMEntity_multiIndex &dofs,ublas::vector<int> &nodes_indices) {
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

PetscErrorCode ForcesAndSurcesCore::getRowNodesIndices(const string &field_name) {
    PetscFunctionBegin;
    PetscErrorCode ierr;
    ierr = getNodesIndices(field_name,const_cast<FENumeredDofMoFEMEntity_multiIndex&>(fe_ptr->get_rows_dofs()),rowNodesIndices); CHKERRQ(ierr);
    PetscFunctionReturn(0);
  }

PetscErrorCode ForcesAndSurcesCore::getColNodesIndices(const string &field_name) {
    PetscFunctionBegin;
    PetscErrorCode ierr;
    ierr = getNodesIndices(field_name,const_cast<FENumeredDofMoFEMEntity_multiIndex&>(fe_ptr->get_cols_dofs()),colNodesIndices); CHKERRQ(ierr);
    PetscFunctionReturn(0);
  }

PetscErrorCode ForcesAndSurcesCore::getTypeIndices(const string &field_name,FENumeredDofMoFEMEntity_multiIndex &dofs,EntityType type,int side_number,ublas::vector<int> &indices) {
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

PetscErrorCode ForcesAndSurcesCore::getTypeIndices(const string &field_name,FENumeredDofMoFEMEntity_multiIndex &dofs,EntityType type,vector<ublas::vector<int> > &indices) {
    PetscFunctionBegin;
    PetscErrorCode ierr;
    SideNumber_multiIndex& side_table = const_cast<SideNumber_multiIndex&>(fe_ptr->get_side_number_table());
    SideNumber_multiIndex::nth_index<2>::type::iterator siit = side_table.get<2>().lower_bound(type);
    SideNumber_multiIndex::nth_index<2>::type::iterator hi_siit = side_table.get<2>().upper_bound(type);
    indices.resize(side_table.get<2>().count(type));
    for(;siit!=hi_siit;siit++) {
      ierr = getTypeIndices(field_name,dofs,type,siit->side_number,indices[siit->side_number]); CHKERRQ(ierr);
    }
    PetscFunctionReturn(0);
  }

PetscErrorCode ForcesAndSurcesCore::getEdgeRowIndices(const string &field_name) {
    PetscFunctionBegin;
    PetscErrorCode ierr;
    ierr = getTypeIndices(field_name,const_cast<FENumeredDofMoFEMEntity_multiIndex&>(fe_ptr->get_rows_dofs()),MBEDGE,rowEdgesIndcies); CHKERRQ(ierr);
    PetscFunctionReturn(0);
  }

PetscErrorCode ForcesAndSurcesCore::getEdgeColIndices(const string &field_name) {
    PetscFunctionBegin;
    PetscErrorCode ierr;
    ierr = getTypeIndices(field_name,const_cast<FENumeredDofMoFEMEntity_multiIndex&>(fe_ptr->get_cols_dofs()),MBEDGE,colEdgesIndcies); CHKERRQ(ierr);
    PetscFunctionReturn(0);
  }

PetscErrorCode ForcesAndSurcesCore::getFacesRowIndices(const string &field_name) {
    PetscFunctionBegin;
    PetscErrorCode ierr;
    ierr = getTypeIndices(field_name,const_cast<FENumeredDofMoFEMEntity_multiIndex&>(fe_ptr->get_rows_dofs()),MBTRI,rowFacesIndcies); CHKERRQ(ierr);
    PetscFunctionReturn(0);
  }

PetscErrorCode ForcesAndSurcesCore::getFacesColIndices(const string &field_name) {
    PetscFunctionBegin;
    PetscErrorCode ierr;
    ierr = getTypeIndices(field_name,const_cast<FENumeredDofMoFEMEntity_multiIndex&>(fe_ptr->get_cols_dofs()),MBTRI,colFacesIndcies); CHKERRQ(ierr);
    PetscFunctionReturn(0);
  }

PetscErrorCode ForcesAndSurcesCore::getTetRowIndices(const string &field_name) {
    PetscFunctionBegin;
    PetscErrorCode ierr;
    ierr = getTypeIndices(field_name,const_cast<FENumeredDofMoFEMEntity_multiIndex&>(fe_ptr->get_rows_dofs()),MBTET,0,rowTetIndcies); CHKERRQ(ierr);
    PetscFunctionReturn(0);
  }

PetscErrorCode ForcesAndSurcesCore::getTetColIndices(const string &field_name) {
    PetscFunctionBegin;
    PetscErrorCode ierr;
    ierr = getTypeIndices(field_name,const_cast<FENumeredDofMoFEMEntity_multiIndex&>(fe_ptr->get_cols_dofs()),MBTET,0,colTetIndcies); CHKERRQ(ierr);
    PetscFunctionReturn(0);
  }

PetscErrorCode ForcesAndSurcesCore::getFaceNodes() {
    PetscFunctionBegin;
    //PetscAttachDebugger();
    ErrorCode rval;
    facesNodes.resize(4,3);
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
      for(int nn = 0;nn<3;nn++) facesNodes(side->side_number,nn) = face_conn[nn]; 
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
    const double *G_X,const double *G_Y,const double *G_Z,const int G_DIM) {
    PetscFunctionBegin;

    PetscErrorCode ierr;

    Nnodes_H1.resize(G_DIM,4);
    ierr = ShapeMBTET(&*Nnodes_H1.data().begin(),G_X,G_Y,G_Z,G_DIM); CHKERRQ(ierr);
    diffNnodes_H1.resize(4,3);
    ierr = ShapeDiffMBTET(&*diffNnodes_H1.data().begin()); CHKERRQ(ierr);

    //edges
    if(edgesSense.size()!=6) {
      SETERRQ(PETSC_COMM_SELF,1,"data inconsistency");
    }
    if(edgesOrder.size()!=6) {
      SETERRQ(PETSC_COMM_SELF,1,"data inconsistency");
    }
    Nedges_H1.resize(6);
    double *_H1edgeN_[6];
    for(int ee = 0;ee<6;ee++) {
      Nedges_H1[ee].resize(G_DIM,NBEDGE_H1(edgesOrder[ee]));
      _H1edgeN_[ee] = &*Nedges_H1[ee].data().begin();
    }
    ierr = H1_EdgeShapeFunctions_MBTET(
      &*edgesSense.data().begin(),&*edgesOrder.data().begin(),&*Nnodes_H1.data().begin(),diffNnodes_H1.data().begin(),_H1edgeN_,NULL,G_DIM); CHKERRQ(ierr);

    //faces
    if(facesOrder.size()!=4) {
      SETERRQ(PETSC_COMM_SELF,1,"data inconsistency");
    }
    Nfaces_H1.resize(4);
    double *_H1faceN_[4];
    for(int ff = 0;ff<4;ff++) {
      Nfaces_H1[ff].resize(G_DIM,NBFACE_H1(facesOrder[ff]));
      _H1faceN_[ff] = &*Nfaces_H1[ff].data().begin();
    }
    if(facesNodes.size1() != 4) {
      SETERRQ(PETSC_COMM_SELF,1,"data inconsistency");
    }
    if(facesNodes.size2() != 3) {
      SETERRQ(PETSC_COMM_SELF,1,"data inconsistency");
    }
    ierr = H1_FaceShapeFunctions_MBTET(
      &*facesNodes.data().begin(),&*facesOrder.data().begin(),
      &*Nnodes_H1.data().begin(),diffNnodes_H1.data().begin(),_H1faceN_,NULL,G_DIM); CHKERRQ(ierr);

    //volume
    Nvolume_H1.resize(G_DIM,NBVOLUME_H1(volumeOrder));
    ierr = H1_VolumeShapeFunctions_MBTET(volumeOrder,&*Nnodes_H1.data().begin(),diffNnodes_H1.data().begin(),&*Nvolume_H1.data().begin(),NULL,G_DIM); CHKERRQ(ierr);

    PetscFunctionReturn(0);
  }
  
PetscErrorCode ForcesAndSurcesCore::shapeTETFunctions_L2(const double *G_X,const double *G_Y,const double *G_Z,const int G_DIM) {
    PetscFunctionBegin;

    PetscErrorCode ierr;

    Nnodes_L2.resize(G_DIM,4);
    ierr = ShapeMBTET(&*Nnodes_L2.data().begin(),G_X,G_Y,G_Z,G_DIM); CHKERRQ(ierr);
    diffNnodes_L2.resize(4,3);
    ierr = ShapeDiffMBTET(&*diffNnodes_L2.data().begin()); CHKERRQ(ierr);

    Nvolume_L2.resize(G_DIM,NBVOLUME_L2(volumeOrder));
    ierr = L2_ShapeFunctions_MBTET(volumeOrder,&*Nnodes_L2.data().begin(),&*diffNnodes_L2.data().begin(),&*Nvolume_L2.data().begin(),NULL,G_DIM); CHKERRQ(ierr);
    PetscFunctionReturn(0);
  }


PetscErrorCode ForcesAndSurcesCore::shapeTRIFunctions_H1(
    ApproximationOrder order,const double *G_X,const double *G_Y,const int G_DIM) {
    PetscFunctionBegin;
    
    SETERRQ(PETSC_COMM_SELF,1,"not implemented");

    PetscFunctionReturn(0);
  }


PetscErrorCode mult_H1_H1::calculate_H1_H1_nonsymetric(double *G_W) {
    PetscFunctionBegin;
    if(B1_Nnodes.size1() != B2_Nnodes.size2()) {
      SETERRQ(PETSC_COMM_SELF,1,"data inconsistency");
    }
    int G_DIM = B1_Nnodes.size1();

    H1_H1_nodes.resize(B1_Nnodes.size2(),B2_Nnodes.size2());
    ublas::noalias(H1_H1_nodes) 
	= ublas::zero_matrix<FieldData>( H1_H1_nodes.size1(),H1_H1_nodes.size2() );

    for(int gg = 0;gg<G_DIM;gg++) {

      cblas_dger(
	CblasRowMajor,B1_Nnodes.size2(),B2_Nnodes.size2(),
	G_W[gg],&B1_Nnodes(gg,0),1,&B2_Nnodes(gg,0),1,
	&*H1_H1_nodes.data().begin(),H1_H1_nodes.size2());

    }
    PetscFunctionReturn(0);
  }

PetscErrorCode mult_H1_H1::calculate_H1_H1_symmetric(double *G_W) {
    PetscFunctionBegin;
    if(B1_Nnodes.size1() != B2_Nnodes.size2()) {
      SETERRQ(PETSC_COMM_SELF,1,"data inconsistency");
    }
    int G_DIM = B1_Nnodes.size1();

    H1_H1_nodes.resize(B1_Nnodes.size2(),B2_Nnodes.size2());
    ublas::noalias(H1_H1_nodes) 
	= ublas::zero_matrix<FieldData>( H1_H1_nodes.size1(),H1_H1_nodes.size2() );

    for(int gg = 0;gg<G_DIM;gg++) {

      cblas_dspr(CblasRowMajor,CblasUpper,
	B1_Nnodes.size2(),G_W[gg],&B1_Nnodes(gg,0),1,
	&*H1_H1_nodes.data().begin());

    }
    PetscFunctionReturn(0);
  }


