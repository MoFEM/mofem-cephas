/** \file moabFEMethod_LowLevelStudent.cpp
 * \brief Myltindex containes, data structures and other low-level functions 
 * 
 * Copyright (C) 2013, Lukasz Kaczmarczyk (likask AT wp.pl) <br>
 *
 * The MoFEM package is copyrighted by Lukasz Kaczmarczyk. 
 * It can be freely used for educational and research purposes 
 * by other institutions. If you use this softwre pleas cite my work. 
 *
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
 * License along with MoFEM. If not, see <http://www.gnu.org/licenses/>
*/

#include<moabFEMethod_LowLevelStudent.hpp>
#include<FEM.h>

namespace MoFEM {

const int debug = 1;

FEMethod_LowLevelStudent::FEMethod_LowLevelStudent(Interface& _moab,int _verbose): FEMethod(_moab),ParentMethod(NULL),verbose(_verbose),fe_ent_ptr(NULL) {
  ShapeMBTET(NTET,G_TET_X1,G_TET_Y1,G_TET_Z1,1);
  ShapeDiffMBTET(diffNTET);
  ShapeMBTRI_GAUSS(NTRI,G_TRI_X1,G_TRI_Y1,1);
  ShapeDiffMBTRI(diffNTRI);
}
FEMethod_LowLevelStudent::~FEMethod_LowLevelStudent() {
  if(ParentMethod!=NULL) {
    delete ParentMethod;
  }
}
PetscErrorCode FEMethod_LowLevelStudent::GetDataView(const string &field_name,EntityType type,int side_number_low,int side_number_hi,FEDofMoFEMEntity_multiIndex_view &dof_view) {
  PetscFunctionBegin;
  dofs_by_Composite::iterator miit = data_multiIndex->get<Composite_mi_tag>().lower_bound(boost::make_tuple(field_name,type,side_number_low));
  dofs_by_Composite::iterator hi_miit = data_multiIndex->get<Composite_mi_tag>().upper_bound(boost::make_tuple(field_name,type,side_number_hi));
  for(;miit!=hi_miit;miit++) {
    dof_view.insert(&*miit);
  }
  PetscFunctionReturn(0);
}
PetscErrorCode FEMethod_LowLevelStudent::GetRowView(const string &field_name,EntityType type,int side_number_low,int side_number_hi,FENumeredDofMoFEMEntity_multiIndex_view &dof_view) {
  PetscFunctionBegin;
  numered_dofs_by_Composite::iterator miit = row_multiIndex->get<Composite_mi_tag>().lower_bound(boost::make_tuple(field_name,type,side_number_low));
  numered_dofs_by_Composite::iterator hi_miit = row_multiIndex->get<Composite_mi_tag>().upper_bound(boost::make_tuple(field_name,type,side_number_hi));
  for(;miit!=hi_miit;miit++) {
    dof_view.insert(&*miit);
  }
  PetscFunctionReturn(0);
}
PetscErrorCode FEMethod_LowLevelStudent::GetColView(const string &field_name,EntityType type,int side_number_low,int side_number_hi,FENumeredDofMoFEMEntity_multiIndex_view &dof_view) {
  PetscFunctionBegin;
  numered_dofs_by_Composite::iterator miit = col_multiIndex->get<Composite_mi_tag>().lower_bound(boost::make_tuple(field_name,type,side_number_low));
  numered_dofs_by_Composite::iterator hi_miit = col_multiIndex->get<Composite_mi_tag>().upper_bound(boost::make_tuple(field_name,type,side_number_hi));
  for(;miit!=hi_miit;miit++) {
    dof_view.insert(&*miit);
  }
  PetscFunctionReturn(0);
}
PetscErrorCode FEMethod_LowLevelStudent::preProcess() {
  PetscFunctionBegin;
  PetscFunctionReturn(0);
}
PetscErrorCode FEMethod_LowLevelStudent::operator()() {
  PetscFunctionBegin;
  PetscFunctionReturn(0);
}
PetscErrorCode FEMethod_LowLevelStudent::postProcess() {
  PetscFunctionBegin;
  PetscFunctionReturn(0);
}

//
template <typename T>
struct UnaryOP_PetscGlobalIdx {
  DofIdx operator()(const T *it) { return FEMethod_LowLevelStudent::UnaryFunction_PetscGlobalIdx(it); }
};
template <typename T>
struct UnaryOP_FieldData {
  FieldData operator()(const T *it) { return FEMethod_LowLevelStudent::UnaryFunction_FieldData(it); }
};
template <typename T>
struct UnaryOP_ApproxRank {
  DofIdx operator()(const T *it) { return FEMethod_LowLevelStudent::UnaryFunction_ApproxRank(it); }
};
template <typename T>
struct UnaryOP_ApproxOrder {
  DofIdx operator()(const T *it) { return FEMethod_LowLevelStudent::UnaryFunction_ApproxOrder(it); }
};
template <typename T>
struct UnaryOP_EntDofIdx {
  DofIdx operator()(const T *it) { return FEMethod_LowLevelStudent::UnaryFunction_EntDofIdx(it); }
};
//
template <typename It,typename M1,typename M2,typename UnaryOp>
PetscErrorCode MapDataTET(It &it,
  M1 &nodes,M2 &edges,M2 &faces,M2 &volume) {
  PetscFunctionBegin;
  ApproximationOrder max_order = it->get_max_order();
  ApproximationRank max_rank = it->get_max_rank();
  int side_number = it->side_number_ptr->side_number;
  int nb_dofs_for_order;
  const MoFEMField* field_ptr = it->get_MoFEMField_ptr();
  const MoFEMEntity* ent_ptr = it->get_MoFEMEntity_ptr();
  switch (it->get_ent_type()) {
    case MBVERTEX:
      assert(side_number>=0);
      assert(side_number<4);
      nb_dofs_for_order = 1;
      nodes[field_ptr].resize(max_rank*4,-1);
      (nodes[field_ptr])[side_number*max_rank + it->get_dof_rank()] =  UnaryOp()(&*it);
      break;
    case MBEDGE:
      assert(side_number>=0);
      assert(side_number<6);
      nb_dofs_for_order = it->forder_edge(max_order);
      edges[ent_ptr].resize(max_rank*nb_dofs_for_order,-1);
      edges[ent_ptr][it->get_EntDofIdx()] = UnaryOp()(&*it);
      break;
    case MBTRI:
      assert(side_number>=0);
      assert(side_number<4);
      nb_dofs_for_order = it->forder_face(max_order);
      faces[ent_ptr].resize(max_rank*nb_dofs_for_order,-1);
      faces[ent_ptr][it->get_EntDofIdx()] = UnaryOp()(&*it);
      break;
    case MBTET:
      nb_dofs_for_order = it->forder_elem(max_order);
      volume[ent_ptr].resize(max_rank*nb_dofs_for_order,-1);
      volume[ent_ptr][it->get_EntDofIdx()] = UnaryOp()(&*it);
      break;
    default:
      SETERRQ(PETSC_COMM_SELF,1,"data inconsistency");
  }
  PetscFunctionReturn(0);
}
//
template <typename It,typename M1,typename M2,typename UnaryOp>
PetscErrorCode MapDataPRISM(It &it,
  M1 &nodes,M2 &edges,M2 &faces,M2 &volume) {
  PetscFunctionBegin;
  ApproximationOrder max_order = it->get_max_order();
  ApproximationRank max_rank = it->get_max_rank();
  int side_number = it->side_number_ptr->side_number;
  int nb_dofs_for_order;
  const MoFEMField* field_ptr = it->get_MoFEMField_ptr();
  const MoFEMEntity* ent_ptr = it->get_MoFEMEntity_ptr();
  switch (it->get_ent_type()) {
    case MBVERTEX:
      assert(side_number>=0);
      assert(side_number<6);
      nb_dofs_for_order = 1;
      nodes[field_ptr].resize(max_rank*6,-1);
      (nodes[field_ptr])[side_number*max_rank + it->get_dof_rank()] =  UnaryOp()(&*it);
      break;
    case MBEDGE:
      assert(side_number>=0);
      assert(side_number<9);
      nb_dofs_for_order = it->forder_edge(max_order);
      edges[ent_ptr].resize(max_rank*nb_dofs_for_order,-1);
      edges[ent_ptr][it->get_EntDofIdx()] = UnaryOp()(&*it);
      break;
    case MBTRI:
      assert(side_number==3||side_number==4);
      nb_dofs_for_order = it->forder_face(max_order);
      faces[ent_ptr].resize(max_rank*nb_dofs_for_order,-1);
      faces[ent_ptr][it->get_EntDofIdx()] = UnaryOp()(&*it);
      break;
    case MBPRISM:
      nb_dofs_for_order = it->forder_elem(max_order);
      volume[ent_ptr].resize(max_rank*nb_dofs_for_order,-1);
      volume[ent_ptr][it->get_EntDofIdx()] = UnaryOp()(&*it);
      break;
    default:
      SETERRQ(PETSC_COMM_SELF,1,"data inconsistency");
  }
  PetscFunctionReturn(0);
}
//
template <typename T>
PetscErrorCode SetMaxOrder(const T &miit,
  vector<int> *e,vector<int> *f,int *v) {
  PetscFunctionBegin;
  switch(miit->get_ent_type()) {
    case MBVERTEX:
    break;
    case MBEDGE: {
      if(e == NULL) SETERRQ(PETSC_COMM_SELF,1,"data inconsistency");
      assert(e->size() > (unsigned int)miit->side_number_ptr->side_number);
      int &a = (*e)[miit->side_number_ptr->side_number];
      a = a < miit->get_max_order() ? miit->get_max_order() : a;
    }
    break;
    case MBTRI: {
      if(f == NULL) SETERRQ(PETSC_COMM_SELF,1,"data inconsistency");
      assert(f->size() > (unsigned int)miit->side_number_ptr->side_number);
      int &a = (*f)[miit->side_number_ptr->side_number];
      a = a < miit->get_max_order() ? miit->get_max_order() : a;
    }
    break;
    case MBTET: {
      if(v == NULL) SETERRQ(PETSC_COMM_SELF,1,"data inconsistency");
      int &a = (*v);
      a = a < miit->get_max_order() ? miit->get_max_order() : a;
    }
    break;
    default:
      SETERRQ(PETSC_COMM_SELF,1,"data inconsistency");
  }
  PetscFunctionReturn(0);
}
PetscErrorCode FEMethod_LowLevelStudent::InitDataStructures() {
  PetscFunctionBegin;
  // node
  row_nodesGlobIndices.clear();
  col_nodesGlobIndices.clear();
  // edge
  row_edgesGlobIndices.clear();
  col_edgesGlobIndices.clear();
  // face
  row_facesGlobIndices.clear();
  col_facesGlobIndices.clear();
  // vol
  row_elemGlobIndices.clear();
  col_elemGlobIndices.clear();
  // row
  row_N_Matrix_nodes.clear();
  row_N_Matrix_edges.clear();
  row_N_Matrix_faces.clear();
  row_N_Matrix_elem.clear();
  //
  isH1 = isHdiv = isHcurl = isL2 = false;
  //
  data_nodes.clear();
  data_edges.clear();
  data_faces.clear();
  data_elem.clear();
  data_at_gauss_pt.clear();
  diff_data_at_gauss_pt.clear();
  //
  fill(maxOrderEdgeH1.begin(),maxOrderEdgeH1.end(),0);
  fill(maxOrderEdgeHdiv.begin(),maxOrderEdgeHdiv.end(),0);
  fill(maxOrderFaceH1.begin(),maxOrderFaceH1.end(),0);
  fill(maxOrderFaceHdiv.begin(),maxOrderFaceHdiv.end(),0);
  fill(maxOrderFaceHcurl.begin(),maxOrderFaceHcurl.end(),0);
  maxOrderElemH1 = maxOrderElemHdiv = maxOrderElemHcurl = maxOrderElemL2 = 0;
  //
  H1edgeN_TRI.clear();
  diffH1edgeN_TRI.clear();
  H1faceN_TRI.clear();
  diffH1faceN_TRI.clear();
  //
  EntityHandle fe_handle = fe_ent_ptr->get_ent();
  // geometry
  int num_nodes;
  rval = moab.get_connectivity(fe_handle,conn,num_nodes,true); CHKERR_PETSC(rval);
  coords.resize(num_nodes*3);
  rval = moab.get_coords(conn,num_nodes,&coords[0]); CHKERR_PETSC(rval);
  switch(fe_ent_ptr->get_ent_type()) {
    case MBTET: 
      // edge
      maxOrderEdgeH1.resize(6);
      maxOrderEdgeHdiv.resize(6);
      // face
      maxOrderFaceH1.resize(4);
      maxOrderFaceHdiv.resize(4);
      maxOrderFaceHcurl.resize(4);
      break;
    case MBPRISM:
      // edge
      maxOrderEdgeH1.resize(9);
      maxOrderEdgeHdiv.resize(9);
      // face
      maxOrderFaceH1.resize(5);
      maxOrderFaceHdiv.resize(5);
      maxOrderFaceHcurl.resize(5);
      break;
    default:
      SETERRQ(PETSC_COMM_SELF,1,"not implemented");
  }
  PetscFunctionReturn(0);
}
PetscErrorCode FEMethod_LowLevelStudent::GlobIndices() {
  PetscFunctionBegin;
  if(fe_ptr==NULL) SETERRQ(PETSC_COMM_SELF,1,"data inconsistency");
  FENumeredDofMoFEMEntity_multiIndex &rows_dofs = const_cast<FENumeredDofMoFEMEntity_multiIndex&>(*row_multiIndex);
  FENumeredDofMoFEMEntity_multiIndex &cols_dofs = const_cast<FENumeredDofMoFEMEntity_multiIndex&>(*col_multiIndex);
  //EntityHandle ent = fe_ptr->get_ent();
  switch (fe_ptr->get_ent_type()) {
    case MBTET: {
	FENumeredDofMoFEMEntity_multiIndex::iterator miit = rows_dofs.begin();
	for(;miit!=rows_dofs.end();miit++) {
	  switch(miit->get_space()) {
	    case H1: {
	      isH1 = true;
	      ierr = SetMaxOrder(miit, &(maxOrderEdgeH1), &(maxOrderFaceH1), &(maxOrderElemH1) ); CHKERRQ(ierr);
	    }
	    break;
	    case Hdiv: {
	      isHdiv = true;
	      ierr = SetMaxOrder(miit, &(maxOrderEdgeHdiv), &(maxOrderFaceHdiv), &(maxOrderElemHdiv) ); CHKERRQ(ierr);
	    }
	    break;
	    case Hcurl: {
	      isHcurl = true;
	      ierr = SetMaxOrder(miit, NULL, &(maxOrderFaceHcurl), &(maxOrderElemHcurl) ); CHKERRQ(ierr);
	    }
	    break;
	    case L2: {
	      ierr = SetMaxOrder(miit, NULL, NULL, &(maxOrderElemL2) ); CHKERRQ(ierr);
	      isL2 = true;
	    }
	    break;
	    default:
	    SETERRQ(PETSC_COMM_SELF,1,"not implemented");
	  }
	  //PetscGlobIndices
	  ierr = MapDataTET<
	    FENumeredDofMoFEMEntity_multiIndex::iterator, 
	    GlobIndices_Type,GlobIndices_EntType,
	    UnaryOP_PetscGlobalIdx<FENumeredDofMoFEMEntity> >(
	      miit, row_nodesGlobIndices, row_edgesGlobIndices,
	      row_facesGlobIndices, row_elemGlobIndices); CHKERRQ(ierr);
	}
	miit = cols_dofs.begin();
	for(;miit!=cols_dofs.end();miit++) {
	  switch(miit->get_space()) {
	    case H1: {
	      isH1 = true;
	      ierr = SetMaxOrder(miit, &(maxOrderEdgeH1), &(maxOrderFaceH1), &(maxOrderElemH1) ); CHKERRQ(ierr);
	    }
	    break;
	    case Hdiv: {
	      isHdiv = true;
	      ierr = SetMaxOrder(miit, &(maxOrderEdgeHdiv), &(maxOrderFaceHdiv), &(maxOrderElemHdiv) ); CHKERRQ(ierr);
	    }
	    break;
	    case Hcurl: {
	      isHcurl = true;
	      ierr = SetMaxOrder(miit, NULL, &(maxOrderFaceHcurl), &(maxOrderElemHcurl) ); CHKERRQ(ierr);
	    }
	    break;
	    case L2: {
	      ierr = SetMaxOrder(miit, NULL, NULL, &(maxOrderElemL2) ); CHKERRQ(ierr);
	      isL2 = true;
	    }
	    break;
	    default:
	    SETERRQ(PETSC_COMM_SELF,1,"not implemented");
	  }
	  //PetscGlobIndices
	  ierr = MapDataTET<
	    FENumeredDofMoFEMEntity_multiIndex::iterator, 
	    GlobIndices_Type,GlobIndices_EntType,
	    UnaryOP_PetscGlobalIdx<FENumeredDofMoFEMEntity> >(
	      miit, col_nodesGlobIndices, col_edgesGlobIndices,
	      col_facesGlobIndices, col_elemGlobIndices); CHKERRQ(ierr);
	}
      }
      break;
    case MBPRISM: {
	FENumeredDofMoFEMEntity_multiIndex::iterator miit = rows_dofs.begin();
	for(;miit!=rows_dofs.end();miit++) {
	  switch(miit->get_space()) {
	    case H1: {
	      isH1 = true;
	      ierr = SetMaxOrder(miit, &(maxOrderEdgeH1), &(maxOrderFaceH1), &(maxOrderElemH1) ); CHKERRQ(ierr);
	    }
	    break;
	    default:
	    SETERRQ(PETSC_COMM_SELF,1,"not implemented");
	  }
	  MapDataPRISM<
	    FENumeredDofMoFEMEntity_multiIndex::iterator, 
	    GlobIndices_Type,GlobIndices_EntType,
	    UnaryOP_PetscGlobalIdx<FENumeredDofMoFEMEntity> >(
	      miit, row_nodesGlobIndices, row_edgesGlobIndices,
	      row_facesGlobIndices, row_elemGlobIndices); CHKERRQ(ierr);
	}
	miit = cols_dofs.begin();
	for(;miit!=cols_dofs.end();miit++) {
	  switch(miit->get_space()) {
	    case H1: {
	      isH1 = true;
	      ierr = SetMaxOrder(miit, &(maxOrderEdgeH1), &(maxOrderFaceH1), &(maxOrderElemH1) ); CHKERRQ(ierr);
	    }
	    break;
	    default:
	    SETERRQ(PETSC_COMM_SELF,1,"not implemented");
	  }
	  MapDataPRISM<
	    FENumeredDofMoFEMEntity_multiIndex::iterator, 
	    GlobIndices_Type,GlobIndices_EntType,
	    UnaryOP_PetscGlobalIdx<FENumeredDofMoFEMEntity> >(
	      miit, col_nodesGlobIndices, col_edgesGlobIndices,
	      col_facesGlobIndices, col_elemGlobIndices); CHKERRQ(ierr);
	}
    }
    break;
    default:
      SETERRQ(PETSC_COMM_SELF,1,"not implemented");
  }
  PetscFunctionReturn(0);
}
PetscErrorCode FEMethod_LowLevelStudent::ParentData(const string &fe_name) {
  PetscFunctionBegin;
  if(ParentMethod == NULL) {
    ParentMethod = new FEMethod_LowLevelStudent(moab,verbose);
    ierr = ParentMethod->set_problem(problem_ptr); CHKERRQ(ierr);
    ierr = ParentMethod->set_dofs_multiIndex(dofs_moabfield); CHKERRQ(ierr);
    ierr = ParentMethod->set_fes_multiIndex(finite_elements); CHKERRQ(ierr);
    ierr = ParentMethod->set_fes_data_multiIndex(finite_elements_data); CHKERRQ(ierr);
    ierr = ParentMethod->set_adjacencies(fem_adjacencies); CHKERRQ(ierr);
    ierr = ParentMethod->preProcess(); CHKERRQ(ierr);
  }
  EntityHandle parent = fe_ptr->get_parent_ent();
  if(verbose>2) {
    PetscPrintf(PETSC_COMM_WORLD,"Parent ent %u\n",parent);
  }
  EntMoFEMFE_multiIndex::index<Composite_mi_tag>::type::iterator 
    miit =  finite_elements_data->get<Composite_mi_tag>().lower_bound(boost::make_tuple(parent,fe_name));
  EntMoFEMFE_multiIndex::index<Composite_mi_tag>::type::iterator 
    hi_miit = finite_elements_data->get<Composite_mi_tag>().upper_bound(boost::make_tuple(parent,fe_name));
  if(distance(miit,hi_miit) > 1) SETERRQ(PETSC_COMM_SELF,1,"data inconsitency");
  for(;miit!=hi_miit;miit++) {
    if(verbose>2) {
      ostringstream ss;
      ss << *miit << endl;
      PetscPrintf(PETSC_COMM_WORLD,ss.str().c_str());
    }
    assert(ParentMethod!=NULL);
    assert(ParentMethod->fe_ptr==NULL);
    ParentMethod->fe_ent_ptr = &*miit;
    ierr = ParentMethod->set_data_multIndex( const_cast<FEDofMoFEMEntity_multiIndex*>(&(miit->data_dofs)) ); CHKERRQ(ierr);
    ierr = ParentMethod->InitDataStructures();
  }
  PetscFunctionReturn(0);
}
PetscErrorCode FEMethod_LowLevelStudent::DataOp() {
  PetscFunctionBegin;
  FEDofMoFEMEntity_multiIndex &data_dofs = const_cast<FEDofMoFEMEntity_multiIndex&>(*data_multiIndex);
  switch (fe_ent_ptr->get_ent_type()) {
    case MBTET: {
	FEDofMoFEMEntity_multiIndex::iterator miit2 = data_dofs.begin();
	if(miit2==data_dofs.end()) SETERRQ1(PETSC_COMM_SELF,1,"Aaaaaaa .... data_dofs size = %u",data_dofs.size());
	for(;miit2!=data_dofs.end();miit2++) {
	  switch(miit2->get_space()) {
	    case H1: {
	      isH1 = true;
	      ierr = SetMaxOrder(miit2, &(maxOrderEdgeH1), &(maxOrderFaceH1), &(maxOrderElemH1) ); CHKERRQ(ierr);
	    }
	    break;
	    case Hdiv: {
	      isHdiv = true;
	      ierr = SetMaxOrder(miit2, &(maxOrderEdgeHdiv), &(maxOrderFaceHdiv), &(maxOrderElemHdiv) ); CHKERRQ(ierr);
	    }
	    break;
	    case Hcurl: {
	      isHcurl = true;
	      ierr = SetMaxOrder(miit2, NULL, &(maxOrderFaceHcurl), &(maxOrderElemHcurl) ); CHKERRQ(ierr);
	    }
	    break;
	    case L2: {
	      ierr = SetMaxOrder(miit2, NULL, NULL, &(maxOrderElemL2) ); CHKERRQ(ierr);
	      isL2 = true;
	    }
	    break;
	    default:
	    SETERRQ(PETSC_COMM_SELF,1,"not implemented");
	  }
	  //Data
	  ierr = MapDataTET<
	    FEDofMoFEMEntity_multiIndex::iterator, 
	    Data_Type,Data_EntType,
	    UnaryOP_FieldData<FEDofMoFEMEntity> >(
	      miit2, data_nodes, data_edges, data_faces, data_elem); CHKERRQ(ierr);
	}
      }
      break;
    case MBPRISM: {
	FEDofMoFEMEntity_multiIndex::iterator miit2 = data_dofs.begin();
	if(miit2==data_dofs.end()) SETERRQ1(PETSC_COMM_SELF,1,"Aaaaaaa .... data_dofs size = %u",data_dofs.size());
	for(;miit2!=data_dofs.end();miit2++) {
	  switch(miit2->get_space()) {
	    case H1: {
	      isH1 = true;
	      ierr = SetMaxOrder(miit2, &(maxOrderEdgeH1), &(maxOrderFaceH1), &(maxOrderElemH1) ); CHKERRQ(ierr);
	    }
	    break;
	    default:
	    SETERRQ(PETSC_COMM_SELF,1,"not implemented");
	  }
	  //Data
	  ierr = MapDataPRISM<
	    FEDofMoFEMEntity_multiIndex::iterator, 
	    Data_Type,Data_EntType,
	    UnaryOP_FieldData<FEDofMoFEMEntity> >(
	      miit2, data_nodes, data_edges, data_faces, data_elem); CHKERRQ(ierr);
	}
      }
      break;
    default:
      SETERRQ(PETSC_COMM_SELF,1,"not implemented");
  }
  PetscFunctionReturn(0);
}
PetscErrorCode FEMethod_LowLevelStudent::ShapeFunctions_TET(vector<double>& _gNTET) {
  PetscFunctionBegin;
  const int cannonical_face_sense_p1[4][3] = { {0,1,3}, {1,2,3}, {0,3,2}/**/, {0,2,1}/**/ }; //secon index is offset (positive sense)
  const int cannonical_face_sense_m1[4][3] = { {0,3,1}, {1,3,2}, {0,2,3}, {0,1,2} }; //second index is offset (negative sense)
  gNTET = _gNTET;
  switch (fe_ent_ptr->get_ent_type()) {
    case MBTET: 
      {
	double invJac[9];
	ierr = ShapeJacMBTET(diffNTET,&*coords.begin(),invJac); CHKERRQ(ierr);
	ierr = Shape_invJac(invJac); CHKERRQ(ierr);
	ierr = ShapeDiffMBTETinvJ(diffNTET,invJac,diffNTETinvJac); CHKERRQ(ierr);
	//
	SideNumber_multiIndex& side_table = const_cast<SideNumber_multiIndex&>(fe_ent_ptr->get_side_number_table());
	// edge
	int _sense_edges_[6];
	if(isH1 || isHdiv) {
	  SideNumber_multiIndex::nth_index<1>::type::iterator siit = side_table.get<1>().lower_bound(boost::make_tuple(MBEDGE,0));
	  SideNumber_multiIndex::nth_index<1>::type::iterator hi_siit = side_table.get<1>().upper_bound(boost::make_tuple(MBEDGE,6));
	  for(;siit!=hi_siit;siit++) {
	    _sense_edges_[siit->side_number] = siit->sense;
	  }
	}
	// face
	int _faces_nodes_[4*3];
	if(isH1 || isHdiv || isHcurl) {
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
	    for(int nn = 0;nn<3;nn++) _faces_nodes_[side->side_number*3+nn] = face_conn[nn]; 
  	    if(debug) {
	      const EntityHandle *conn_face;
	      int num_nodes;
	      rval = moab.get_connectivity(side->ent,conn_face,num_nodes,true); CHKERR_PETSC(rval);
	      assert(num_nodes==3);
	      assert(conn_face[0] == conn[face_conn[0]]);
	      assert(conn_face[1] == conn[face_conn[1]]);
	      assert(conn_face[2] == conn[face_conn[2]]);
	    }
	  }
	}
	gNTET_dim = gNTET.size()/4;
	if(isH1) {
	  int _edges_order_[6];
	  copy(maxOrderEdgeH1.begin(),maxOrderEdgeH1.end(),&_edges_order_[0]);
	  int _faces_order_[4];
	  copy(maxOrderFaceH1.begin(),maxOrderFaceH1.end(),&_faces_order_[0]);
	  // edge H1
	  H1edgeN.resize(6);
	  diffH1edgeN.resize(6);
	  diffH1edgeNinvJac.resize(6);
	  double *_H1edgeN_[6],*_diffH1edgeN_[6],*_diffH1edgeNinvJac_[6];
	  for(int ee = 0;ee<6;ee++) {
	    H1edgeN[ee].resize(gNTET_dim*NBEDGE_H1(max_ApproximationOrder));
	    diffH1edgeN[ee].resize(3*gNTET_dim*NBEDGE_H1(max_ApproximationOrder));
	    diffH1edgeNinvJac[ee].resize(3*gNTET_dim*NBEDGE_H1(max_ApproximationOrder));
	    _H1edgeN_[ee] = &(H1edgeN[ee])[0];
	    _diffH1edgeN_[ee] = &(diffH1edgeN[ee])[0];
	    _diffH1edgeNinvJac_[ee] = &(diffH1edgeNinvJac[ee])[0]; 
	  }
	  ierr = H1_EdgeShapeFunctions_MBTET(_sense_edges_,_edges_order_,&gNTET[0],diffNTET,_H1edgeN_,_diffH1edgeN_,gNTET_dim); CHKERRQ(ierr);
	  ierr = H1_EdgeShapeDiffMBTETinvJ(_edges_order_,_edges_order_,_diffH1edgeN_,invJac,_diffH1edgeNinvJac_,gNTET_dim); CHKERRQ(ierr);
	  //copy(diffH1edgeNinvJac[0].begin(),diffH1edgeNinvJac[0].end(),ostream_iterator<double>(cerr," ")); cerr << endl;
	  // face H1
	  H1faceN.resize(4);
	  diffH1faceN.resize(4);
	  diffH1faceNinvJac.resize(4);
	  double *_H1faceN_[4],*_diffH1faceN_[4],*_diffH1faceNinvJac_[4];
	  for(int ff = 0;ff<4;ff++) {
	    H1faceN[ff].resize(gNTET_dim*NBFACE_H1(max_ApproximationOrder));
	    diffH1faceN[ff].resize(3*gNTET_dim*NBFACE_H1(max_ApproximationOrder));
	    diffH1faceNinvJac[ff].resize(3*gNTET_dim*NBFACE_H1(max_ApproximationOrder));
	    _H1faceN_[ff] = &(H1faceN[ff])[0];
	    _diffH1faceN_[ff] = &(diffH1faceN[ff])[0];
	    _diffH1faceNinvJac_[ff] = &(diffH1faceNinvJac[ff])[0];
	  }
	  ierr = H1_FaceShapeFunctions_MBTET(_faces_nodes_,_faces_order_,&gNTET[0],diffNTET,_H1faceN_,_diffH1faceN_,gNTET_dim); CHKERRQ(ierr);
	  ierr = H1_FaceShapeDiffMBTETinvJ(_faces_order_,_faces_order_,_diffH1faceN_,invJac,_diffH1faceNinvJac_,gNTET_dim); CHKERRQ(ierr);
	  // vol H1
	  H1elemN.resize(gNTET_dim*NBVOLUME_H1(max_ApproximationOrder));
	  diffH1elemN.resize(3*gNTET_dim*NBVOLUME_H1(max_ApproximationOrder));
	  diffH1elemNinvJac.resize(3*gNTET_dim*NBVOLUME_H1(max_ApproximationOrder));
	  ierr = H1_VolumeShapeFunctions_MBTET(maxOrderElemH1,&gNTET[0],diffNTET,&H1elemN[0],&diffH1elemN[0],gNTET_dim);  CHKERRQ(ierr);
	  ierr = H1_VolumeShapeDiffMBTETinvJ(maxOrderElemH1,maxOrderElemH1,&diffH1elemN[0],invJac,&diffH1elemNinvJac[0],gNTET_dim); CHKERRQ(ierr);
	}
	if(isL2) {
	  // vol L2
	  L2elemN.resize(gNTET_dim*NBVOLUME_L2(max_ApproximationOrder));
	  diffL2elemN.resize(3*gNTET_dim*NBVOLUME_L2(max_ApproximationOrder));
	  diffL2elemNinvJac.resize(3*gNTET_dim*NBVOLUME_L2(max_ApproximationOrder));
	  ierr = L2_ShapeFunctions_MBTET(maxOrderElemL2,&gNTET[0],diffNTET,&L2elemN[0],&diffL2elemN[0],gNTET_dim); CHKERRQ(ierr);
	  //copy(gNTET.begin(),gNTET.end(),ostream_iterator<double>(cout, " ")); cout << "gNTET" << endl;
	  //copy(&L2elemN[0],&L2elemN[gNTET_dim*NBVOLUME_L2(maxOrderElemL2)],ostream_iterator<double>(cout, " ")); cout << "L2elemN" << endl;
	  ierr = L2_VolumeShapeDiffMBTETinvJ(maxOrderElemL2,maxOrderElemL2,&diffL2elemN[0],invJac,&diffL2elemNinvJac[0],gNTET_dim); CHKERRQ(ierr);
	}
      } 
      break;
      case MBPRISM: {
	SETERRQ(PETSC_COMM_SELF,1,"Aaaa... not implemented yet");
      } 
      break;
     default:
      SETERRQ(PETSC_COMM_SELF,1,"not implemented");
  }
  PetscFunctionReturn(0);
}
PetscErrorCode FEMethod_LowLevelStudent::ShapeFunctions_PRISM(vector<double>& _gNTRI_) {
  PetscFunctionBegin;
  gNTRI = _gNTRI_;
  gNTRI_dim = gNTRI.size()/3;
  switch (fe_ent_ptr->get_ent_type()) {
    case MBPRISM: {
      SideNumber_multiIndex& side_table = const_cast<SideNumber_multiIndex&>(fe_ent_ptr->get_side_number_table());
      SideNumber_multiIndex::nth_index<1>::type::iterator siit3 = side_table.get<1>().find(boost::make_tuple(MBTRI,3));
      SideNumber_multiIndex::nth_index<1>::type::iterator siit4 = side_table.get<1>().find(boost::make_tuple(MBTRI,4));
      if(siit3==side_table.get<1>().end()) SETERRQ(PETSC_COMM_SELF,1,"data inconsistency");
      if(siit4==side_table.get<1>().end()) SETERRQ(PETSC_COMM_SELF,1,"data inconsistency");
      int num_nodes;
      const EntityHandle *conn_face3;
      rval = moab.get_connectivity(siit3->ent,conn_face3,num_nodes,true); CHKERR_PETSC(rval);
      assert(num_nodes==3);
      const EntityHandle *conn_face4;
      rval = moab.get_connectivity(siit4->ent,conn_face4,num_nodes,true); CHKERR_PETSC(rval);
      assert(num_nodes==3);
      if(isH1) {
	//nodes
	gNTRIonPRISM.resize(6*gNTRI_dim);
	try {
	  for(int gg = 0;gg<gNTRI_dim;gg++) {
	    for(int nn = 0;nn<3;nn++) {
	      gNTRIonPRISM[fe_ent_ptr->get_side_number_ptr(moab,conn_face3[nn])->side_number] = +gNTRI[gg*3+nn]; 
	      gNTRIonPRISM[fe_ent_ptr->get_side_number_ptr(moab,conn_face4[nn])->side_number] = -gNTRI[gg*3+nn]; 
	    }
	  }
	} catch (const char* msg) {
	  SETERRQ(PETSC_COMM_SELF,1,msg);
	} 
	//edges
	int _edge_sense3_[3],_edge_sense4_[3];
	EntityHandle edges3[3],edges4[3];
	for(int ee = 0;ee<3;ee++) {
	  rval = moab.side_element(siit3->ent,1,ee,edges3[ee]); CHKERR_PETSC(rval);
	  rval = moab.side_element(siit4->ent,1,ee,edges4[ee]); CHKERR_PETSC(rval);
	  int side_number,offset;
	  rval = moab.side_number(siit3->ent,edges3[ee],side_number,_edge_sense3_[ee],offset); CHKERR_PETSC(rval);
	  assert(side_number == ee);
	  rval = moab.side_number(siit4->ent,edges4[ee],side_number,_edge_sense4_[ee],offset); CHKERR_PETSC(rval);
	  assert(side_number == ee);
	}
	H1edgeN.resize(9);
	diffH1edgeN.resize(0);
	double *_H1edgeN3_[3],*_H1edgeN4_[3];
	int _edge_order3_[3],_edge_order4_[3];
	for(int ee = 0;ee<3;ee++) {
	  H1edgeN[ee].resize(gNTRI_dim*NBEDGE_H1(max_ApproximationOrder));
	  H1edgeN[ee+6].resize(gNTRI_dim*NBEDGE_H1(max_ApproximationOrder));
	  _H1edgeN3_[ee] = &(H1edgeN[ee][0]);
	  _H1edgeN4_[ee] = &(H1edgeN[ee+6][0]);
	  _edge_order3_[ee] = maxOrderEdgeH1[ee];
	  _edge_order4_[ee] = maxOrderEdgeH1[ee+6];
	}
	ierr = H1_EdgeShapeFunctions_MBTRI(_edge_sense3_,_edge_order3_,&gNTRI[0],diffNTRI,_H1edgeN3_,NULL,gNTRI_dim); CHKERRQ(ierr);
	ierr = H1_EdgeShapeFunctions_MBTRI(_edge_sense4_,_edge_order4_,&gNTRI[0],diffNTRI,_H1edgeN4_,NULL,gNTRI_dim); CHKERRQ(ierr);
	for(int ee = 0;ee<3;ee++) {
	  cblas_dscal(NBEDGE_H1(_edge_order4_[ee]),-1,_H1edgeN4_[ee],1);
	}
	//faces
	int _face_order3_ = maxOrderFaceH1[siit3->side_number];
	int _face_order4_ = maxOrderFaceH1[siit4->side_number];
	H1faceN.resize(5);
	diffH1faceN.resize(0);
	H1faceN[3].resize(NBFACE_H1(max_ApproximationOrder)*gNTRI_dim);
	H1faceN[4].resize(NBFACE_H1(max_ApproximationOrder)*gNTRI_dim);
	double *_faceN3_ = &(H1faceN[3][0]);
	double *_faceN4_ = &(H1faceN[4][0]);
	ierr = H1_FaceShapeFunctions_MBTRI(_face_order3_,&gNTRI[0],diffNTRI,_faceN3_,NULL,gNTRI_dim); CHKERRQ(ierr);
	ierr = H1_FaceShapeFunctions_MBTRI(_face_order4_,&gNTRI[0],diffNTRI,_faceN4_,NULL,gNTRI_dim); CHKERRQ(ierr);
	cblas_dscal(NBFACE_H1(_face_order4_),-1,_faceN4_,1);
      }
    }
    break;
    default:
      SETERRQ(PETSC_COMM_SELF,1,"not implemented");
  }
  PetscFunctionReturn(0);
}
PetscErrorCode FEMethod_LowLevelStudent::get_ShapeFunction(
    vector<const double*> *shape_by_gauss_pt,
    vector<const double*> *diff_shape_by_gauss_pt,
    const MoFEMField* field_ptr,EntityType type,int side_number) {
  PetscFunctionBegin;
  switch (fe_ent_ptr->get_ent_type()) {
    case MBTET: {
      int gg = 0;
      if(shape_by_gauss_pt!=NULL) {
	shape_by_gauss_pt->resize(gNTET_dim,NULL);
      }
      if(diff_shape_by_gauss_pt!=NULL) {
	diff_shape_by_gauss_pt->resize(gNTET_dim,NULL);
      }
      switch (field_ptr->get_space()) {
        case H1: {
          switch (type) {
	    case MBVERTEX: {
	      if(side_number != -1) SETERRQ(PETSC_COMM_SELF,1,"data inconsistency");
	      for(;gg<gNTET_dim;gg++) {
		if(shape_by_gauss_pt!=NULL) {
		  (*shape_by_gauss_pt)[gg] = &gNTET[4*gg];
		}
		if(diff_shape_by_gauss_pt!=NULL) {
		  (*diff_shape_by_gauss_pt)[gg] = &diffNTETinvJac[0];
		}
	      }
	    }
	    break;
	    case MBEDGE: {
	      if(side_number < 0) SETERRQ(PETSC_COMM_SELF,1,"data inconsistency");
	      if(side_number > 6) SETERRQ(PETSC_COMM_SELF,1,"data inconsistency");
	      for(;gg<gNTET_dim;gg++) {
		if(shape_by_gauss_pt!=NULL) {
		  (*shape_by_gauss_pt)[gg] = &((H1edgeN[side_number])[gg*NBEDGE_H1(maxOrderEdgeH1[side_number])]);
		}
		if(diff_shape_by_gauss_pt!=NULL) {
		  (*diff_shape_by_gauss_pt)[gg] = &((diffH1edgeNinvJac[side_number])[3*gg*NBEDGE_H1(maxOrderEdgeH1[side_number])]);
		}
	      }
	    }
	    break;
	    case MBTRI: {
	      if(side_number < 0) SETERRQ(PETSC_COMM_SELF,1,"data inconsistency");
	      if(side_number > 4) SETERRQ(PETSC_COMM_SELF,1,"data inconsistency");
	      for(;gg<gNTET_dim;gg++) {
		if(shape_by_gauss_pt!=NULL) {
		  (*shape_by_gauss_pt)[gg] = &((H1faceN[side_number])[gg*NBFACE_H1(maxOrderFaceH1[side_number])]);
		}
  		if(diff_shape_by_gauss_pt!=NULL) {
		  (*diff_shape_by_gauss_pt)[gg] = &((diffH1faceNinvJac[side_number])[3*gg*NBFACE_H1(maxOrderFaceH1[side_number])]);
		}
	      }
	    }
	    break;
	    case MBTET: {
	      if(side_number != -1) SETERRQ(PETSC_COMM_SELF,1,"data inconsistency");
	      for(;gg<gNTET_dim;gg++) {
		if(shape_by_gauss_pt!=NULL) {
		  (*shape_by_gauss_pt)[gg] = &((H1elemN)[gg*NBVOLUME_H1(maxOrderElemH1)]);
		}
  		if(diff_shape_by_gauss_pt!=NULL) {
		  (*diff_shape_by_gauss_pt)[gg] = &(diffH1elemNinvJac[3*gg*NBVOLUME_H1(maxOrderElemH1)]);
		}
	      }
	    }
	    break;  
	    default:
	      SETERRQ(PETSC_COMM_SELF,1,"not implemented");
	  }
	}
	break;
	case L2: {
	  if(type != MBTET) SETERRQ(PETSC_COMM_SELF,1,"data inconsistency");
      	  for(;gg<gNTET_dim;gg++) {
    	    if(shape_by_gauss_pt!=NULL) {
	      (*shape_by_gauss_pt)[gg] = &((L2elemN)[gg*NBVOLUME_L2(maxOrderElemL2)]);
	    }
	    if(diff_shape_by_gauss_pt!=NULL) {
	      (*diff_shape_by_gauss_pt)[gg] = &(diffL2elemNinvJac[3*gg*NBVOLUME_L2(maxOrderElemL2)]);
	    }
	  }
	}
	break;
        default:
          SETERRQ(PETSC_COMM_SELF,1,"not implemented");
      }
      break;
    }
    break;
    case MBPRISM: {
      if(shape_by_gauss_pt!=NULL) {
	shape_by_gauss_pt->resize(gNTRI_dim,NULL);
      }
      if(diff_shape_by_gauss_pt!=NULL) {
	SETERRQ(PETSC_COMM_SELF,1,"Aaaaa... not implemented yet");
      }
      int gg = 0;
      switch (field_ptr->get_space()) {
        case H1: {
          switch (type) {
	    case MBVERTEX: {
	      if(side_number != -1) SETERRQ(PETSC_COMM_SELF,1,"data inconsistency");
	      for(;gg<gNTRI_dim;gg++) {
		if(shape_by_gauss_pt!=NULL) {
		  (*shape_by_gauss_pt)[gg] = &gNTRIonPRISM[6*gg];
		}
	      }
	    }
	    break;
	    case MBEDGE: {
	      if(side_number < 0) SETERRQ(PETSC_COMM_SELF,1,"data inconsistency");
	      if(side_number > 8) SETERRQ(PETSC_COMM_SELF,1,"data inconsistency");
	      if(side_number > 2 && side_number < 6) SETERRQ(PETSC_COMM_SELF,1,"data inconsistency");
	      for(;gg<gNTRI_dim;gg++) {
		if(shape_by_gauss_pt!=NULL) {
		  (*shape_by_gauss_pt)[gg] = &((H1edgeN[side_number])[gg*NBEDGE_H1(maxOrderEdgeH1[side_number])]);
		}
	      }
	    }
	    break;
	    case MBTRI: {
	      if(side_number != 3 && side_number != 4) SETERRQ(PETSC_COMM_SELF,1,"data inconsistency");
	      for(;gg<gNTRI_dim;gg++) {
		if(shape_by_gauss_pt!=NULL) {
		  (*shape_by_gauss_pt)[gg] = &((H1faceN[side_number])[gg*NBFACE_H1(maxOrderFaceH1[side_number])]);
		}
	      }
	    }
	    break;
	    default:
	      SETERRQ(PETSC_COMM_SELF,1,"not implemented");
	  }
	}
	break;
        default:
          SETERRQ(PETSC_COMM_SELF,1,"not implemented");
      }
    }
    break;
    default:
      SETERRQ(PETSC_COMM_SELF,1,"not implemented");
  }
  PetscFunctionReturn(0);
}
PetscErrorCode FEMethod_LowLevelStudent::Data_at_GaussPoints() {
  PetscFunctionBegin;
  unsigned int g_dim,nb_Ns;
  switch (fe_ent_ptr->get_ent_type()) {
    case MBTET:
      g_dim = gNTET.size()/4;
      nb_Ns = 4;
      break;
    case MBPRISM:
      g_dim = gNTRI.size()/3;
      nb_Ns = 6;
      break;
    default:
      SETERRQ(PETSC_COMM_SELF,1,"not implemented yet");
  }
  data_at_gauss_pt.clear();
  // nodes
  for(Data_Type::iterator dit = data_nodes.begin();dit!=data_nodes.end();dit++) {
    const MoFEMField* field_ptr = dit->first;
    const string &field_name = field_ptr->get_name();
    ublas::vector<FieldData> &dof_data = dit->second;
    vector<const double*> shape_by_gauss_pt;
    ierr = get_ShapeFunction(&shape_by_gauss_pt,NULL,field_ptr,MBVERTEX); CHKERRQ(ierr);
    assert(shape_by_gauss_pt.size()==g_dim);
    vector< ublas::vector<FieldData> > &data = data_at_gauss_pt[field_name];
    data.resize(g_dim);
    int rank = field_ptr->get_max_rank();
    unsigned int gg = 0;
    for(;gg<g_dim;gg++) {
      data[gg].resize(rank);
      int rr = 0;
      for(;rr<rank;rr++) {
	if(shape_by_gauss_pt[gg] == NULL) SETERRQ(PETSC_COMM_SELF,1,"data inconsitencies");
	(data[gg])[rr] += cblas_ddot(nb_Ns,shape_by_gauss_pt[gg],1,&dof_data[rr],rank);
      }
    }
  }
  // edges // faces // volumes
  Data_EntType* F[] = { &data_edges, &data_faces, &data_elem };
  for(int ss = 0;ss<3;ss++) {
    for(Data_EntType::iterator dit = F[ss]->begin();dit!=F[ss]->end();dit++) {
      const MoFEMEntity* ent_ptr = dit->first;
      const MoFEMField* field_ptr = ent_ptr->get_MoFEMField_ptr();
      const string &field_name = field_ptr->get_name();
      ublas::vector<FieldData> &dof_data = dit->second;
      vector<const double*> shape_by_gauss_pt;
      if(ss<2) {
	try {
	  int side_number = fe_ent_ptr->get_side_number_ptr(moab,ent_ptr->get_ent())->side_number;
	  ierr = get_ShapeFunction(&shape_by_gauss_pt,NULL,field_ptr,ent_ptr->get_ent_type(),side_number); CHKERRQ(ierr);
	} catch (const char* msg) {
	  SETERRQ(PETSC_COMM_SELF,1,msg);
	}
      } else {
	ierr = get_ShapeFunction(&shape_by_gauss_pt,NULL,field_ptr,ent_ptr->get_ent_type()); CHKERRQ(ierr);
      }
      assert(shape_by_gauss_pt.size()==g_dim);
      vector< ublas::vector<FieldData> > &data = data_at_gauss_pt[field_name];
      data.resize(g_dim);
      unsigned int rank = field_ptr->get_max_rank();
      unsigned int order = ent_ptr->get_max_order();
      unsigned int nb_dofs = ent_ptr->forder(order);
      if(nb_dofs == 0) continue;
      if(dof_data.size()/rank != nb_dofs) SETERRQ(PETSC_COMM_SELF,1,"data inconsitencies");
      if(shape_by_gauss_pt.size()/g_dim > nb_dofs) SETERRQ(PETSC_COMM_SELF,1,"data inconsitencies");
      unsigned int gg = 0;
      for(;gg<g_dim;gg++) {
	data[gg].resize(rank);
        unsigned int rr = 0;
        for(;rr<rank;rr++) {
	  if(shape_by_gauss_pt[gg] == NULL) SETERRQ(PETSC_COMM_SELF,1,"data inconsitencies");
	  (data[gg])[rr] += cblas_ddot(nb_dofs,shape_by_gauss_pt[gg],1,&dof_data[rr],rank);
        }
      }
    }
  }
  PetscFunctionReturn(0);
}
PetscErrorCode FEMethod_LowLevelStudent::DiffData_at_GaussPoints() {
  PetscFunctionBegin;
  unsigned int g_dim,nb_Ns;
  switch (fe_ent_ptr->get_ent_type()) {
    case MBTET:
      g_dim = gNTET.size()/4;
      nb_Ns = 4;
      break;
    default:
      SETERRQ(PETSC_COMM_SELF,1,"not implemented yet");
  }
  diff_data_at_gauss_pt.clear();
  // nodes
  for(Data_Type::iterator dit = data_nodes.begin();dit!=data_nodes.end();dit++) {
    const MoFEMField* field_ptr = dit->first;
    const string &field_name = field_ptr->get_name();
    int dim = 0;
    switch(field_ptr->get_space()) {
      case H1:
      case Hdiv:
      case Hcurl:
      case L2:
	dim = 3;
	break;
      case L2_2D:
	dim = 2;
	break;
      default:
	SETERRQ(PETSC_COMM_SELF,1,"not implemented yet");
    }
    ublas::vector<FieldData> &dof_data = dit->second;
    vector<const double*> diff_shape_by_gauss_pt;
    ierr = get_ShapeFunction(NULL,&diff_shape_by_gauss_pt,field_ptr,MBVERTEX); CHKERRQ(ierr);
    assert(diff_shape_by_gauss_pt.size()==g_dim);
    vector<ublas::matrix<FieldData> > &diff_data = diff_data_at_gauss_pt[field_name];
    diff_data.resize(g_dim);
    int rank = field_ptr->get_max_rank();
    for(unsigned int gg = 0;gg<g_dim;gg++) {
      diff_data[gg] = ublas::zero_matrix<FieldData>(rank,dim);
    }
    unsigned int gg = 0;
    for(;gg<g_dim;gg++) {
      int dd = 0;
      for(;dd<dim;dd++) {
	int rr = 0;
	for(;rr<rank;rr++) {
	  if(diff_shape_by_gauss_pt[gg] == NULL) SETERRQ(PETSC_COMM_SELF,1,"data inconsitencies");
	  (diff_data[gg])(rr,dd) += cblas_ddot(nb_Ns,&(diff_shape_by_gauss_pt[gg])[dd],dim,&dof_data[rr],rank);
	}
      }
    }
  }
  // edges // faces // folumes
  Data_EntType* F[] = { &data_edges, &data_faces, &data_elem };
  for(int ss = 0;ss<3;ss++) {
    for(Data_EntType::iterator dit = F[ss]->begin();dit!=F[ss]->end();dit++) {
      const MoFEMEntity* ent_ptr = dit->first;
      const MoFEMField* field_ptr = ent_ptr->get_MoFEMField_ptr();
      const string &field_name = field_ptr->get_name();
      int dim = 0;
      switch(field_ptr->get_space()) {
	case H1:
	case Hdiv:
	case Hcurl:
	case L2:
	  dim = 3;
	  break;
	case L2_2D:
	  dim = 2;
	  break;
	default:
	  SETERRQ(PETSC_COMM_SELF,1,"not implemented yet");
      }
      ublas::vector<FieldData> &dof_data = dit->second;
      vector<const double*> diff_shape_by_gauss_pt;
      if(ss<2) {
	try {
	  int side_number = fe_ent_ptr->get_side_number_ptr(moab,ent_ptr->get_ent())->side_number;
	  ierr = get_ShapeFunction(NULL,&diff_shape_by_gauss_pt,field_ptr,ent_ptr->get_ent_type(),side_number); CHKERRQ(ierr);
	} catch (const char* msg) {
	  SETERRQ(PETSC_COMM_SELF,1,msg);
	}
      } else {
	ierr = get_ShapeFunction(NULL,&diff_shape_by_gauss_pt,field_ptr,ent_ptr->get_ent_type()); CHKERRQ(ierr);
      }
      assert(diff_shape_by_gauss_pt.size()==g_dim);
      vector<ublas::matrix<FieldData> > &diff_data = diff_data_at_gauss_pt[field_name];
      diff_data.resize(g_dim);
      unsigned int rank = field_ptr->get_max_rank();
      unsigned int order = ent_ptr->get_max_order();
      unsigned int nb_dofs = ent_ptr->forder(order);
      if(nb_dofs == 0) continue;
      if(dof_data.size()/rank != nb_dofs) SETERRQ(PETSC_COMM_SELF,1,"data inconsitencies");
      if(diff_shape_by_gauss_pt.size()/(dim*g_dim) > nb_dofs) SETERRQ(PETSC_COMM_SELF,1,"data inconsitencies");
      unsigned int gg = 0;
      for(;gg<g_dim;gg++) {
	diff_data[gg].resize(rank,dim);
	int dd = 0;
	for(;dd<dim;dd++) {
	  unsigned int rr = 0;
	  for(;rr<rank;rr++) {
	    if(diff_shape_by_gauss_pt[gg] == NULL) SETERRQ(PETSC_COMM_SELF,1,"data inconsitencies");
	    (diff_data[gg])(rr,dd) += cblas_ddot(nb_dofs,&(diff_shape_by_gauss_pt[gg])[dd],dim,&dof_data[rr],rank);
	  }
	}
      }
    }
  }
  PetscFunctionReturn(0);
}
PetscErrorCode FEMethod_LowLevelStudent::GetNMatrix_at_GaussPoint(
    GlobIndices_Type& nodesGlobIndices, GlobIndices_EntType& edgesGlobIndices,
    GlobIndices_EntType& facesGlobIndices, GlobIndices_EntType& volumeGlobIndices,
    N_Matrix_Type& N_Matrix_nodes,
    N_Matrix_EntType& N_Matrix_edges,
    N_Matrix_EntType& N_Matrix_faces,
    N_Matrix_EntType& N_Matrix_elem) {
  PetscFunctionBegin;
  unsigned int g_dim,nb_Ns;
  switch (fe_ent_ptr->get_ent_type()) {
    case MBTET:
      g_dim = gNTET.size()/4;
      nb_Ns = 4;
      break;
    case MBPRISM:
      g_dim = gNTRI.size()/3;
      nb_Ns = 6;
      break;
    default:
      SETERRQ(PETSC_COMM_SELF,1,"not implemented yet");
  }
  // nodes
  GlobIndices_Type::iterator nit = nodesGlobIndices.begin();
  for(;nit!=nodesGlobIndices.end();nit++) {
    const MoFEMField* field_ptr = nit->first;
    vector<const double*> shape_by_gauss_pt;
    ierr = get_ShapeFunction(&shape_by_gauss_pt,NULL,field_ptr,MBVERTEX); CHKERRQ(ierr);
    assert(shape_by_gauss_pt.size()==g_dim);
    vector< ublas::matrix<FieldData> > &data = N_Matrix_nodes[field_ptr];
    data.resize(g_dim);
    int rank = field_ptr->get_max_rank();
    if(rank*nb_Ns!=nit->second.size()) SETERRQ(PETSC_COMM_SELF,1,"data inconsitency");
    unsigned int gg = 0;
    for(;gg<g_dim;gg++) {
      ublas::matrix<FieldData> &mat = data[gg];
      mat.resize(rank,rank*nb_Ns);
      mat = ublas::zero_matrix<FieldData>(rank,rank*nb_Ns);
      int rr = 0;
      for(;rr<rank;rr++) {
	if(shape_by_gauss_pt[gg] == NULL) SETERRQ(PETSC_COMM_SELF,1,"data inconsitency");
	ublas::matrix_row<ublas::matrix<FieldData> > mr(mat,rr);
	for(unsigned int jj = 0;jj<nb_Ns;jj++) mr(rank*jj + rr) = (shape_by_gauss_pt[gg])[jj];
      }
    }
  }
  // edges // faces // volumes
  GlobIndices_EntType* F[] = { &edgesGlobIndices, &facesGlobIndices, &volumeGlobIndices };
  N_Matrix_EntType* FF[] = {  &N_Matrix_edges, &N_Matrix_faces, &N_Matrix_elem };
  for(int ss = 0;ss<3;ss++) {
    for(GlobIndices_EntType::iterator dit = F[ss]->begin();dit!=F[ss]->end();dit++) {
      const MoFEMEntity* ent_ptr = dit->first;
      const MoFEMField* field_ptr = ent_ptr->get_MoFEMField_ptr();
      vector<const double*> shape_by_gauss_pt;
      if(ss<2) {
	try {
	  int side_number = fe_ent_ptr->get_side_number_ptr(moab,ent_ptr->get_ent())->side_number;
	  ierr = get_ShapeFunction(&shape_by_gauss_pt,NULL,field_ptr,ent_ptr->get_ent_type(),side_number); CHKERRQ(ierr);
	} catch (const char* msg) {
	  SETERRQ(PETSC_COMM_SELF,1,msg);
	}
      } else {
	ierr = get_ShapeFunction(&shape_by_gauss_pt,NULL,field_ptr,ent_ptr->get_ent_type()); CHKERRQ(ierr);
      }
      vector<ublas::matrix<FieldData> > &data = (*FF[ss])[ent_ptr];
      data.resize(g_dim);
      int rank = field_ptr->get_max_rank();
      int order = ent_ptr->get_max_order();
      unsigned int nb_dofs = rank*ent_ptr->forder(order); 
      if(nb_dofs!=dit->second.size()) SETERRQ(PETSC_COMM_SELF,1,"data inconsitency");
      unsigned int gg = 0;
      for(;gg<g_dim;gg++) {
	if(shape_by_gauss_pt[gg] == NULL) SETERRQ(PETSC_COMM_SELF,1,"data inconsitency");
	ublas::matrix<FieldData> &mat = data[gg];
	mat.resize(rank,nb_dofs);
	mat = ublas::zero_matrix<FieldData>(rank,nb_dofs);
	for(int rr = 0;rr<rank;rr++) {
	  ublas::matrix_row<ublas::matrix<FieldData> > mr(mat,rr);
	  for(int jj = 0;jj<ent_ptr->forder(order);jj++) {
	    mr(rank*jj + rr) = (shape_by_gauss_pt[gg])[jj];
	  }
	}
      }
    }
  }
  PetscFunctionReturn(0);
}
PetscErrorCode FEMethod_LowLevelStudent::GetRowNMatrix_at_GaussPoint() {
  PetscFunctionBegin;
  ierr = GetNMatrix_at_GaussPoint(
    row_nodesGlobIndices,row_edgesGlobIndices,
    row_facesGlobIndices,row_elemGlobIndices,
    row_N_Matrix_nodes, row_N_Matrix_edges,
    row_N_Matrix_faces, row_N_Matrix_elem); CHKERRQ(ierr);
  PetscFunctionReturn(0);
}
PetscErrorCode FEMethod_LowLevelStudent::GetColNMatrix_at_GaussPoint() {
  PetscFunctionBegin;
  ierr = GetNMatrix_at_GaussPoint(
    col_nodesGlobIndices,col_edgesGlobIndices,
    col_facesGlobIndices,col_elemGlobIndices,
    col_N_Matrix_nodes, col_N_Matrix_edges,
    col_N_Matrix_faces, col_N_Matrix_elem); CHKERRQ(ierr);
  PetscFunctionReturn(0);
}
PetscErrorCode FEMethod_LowLevelStudent::GetDiffNMatrix_at_GaussPoint(
    GlobIndices_Type& nodesGlobIndices, GlobIndices_EntType& edgesGlobIndices,
    GlobIndices_EntType& facesGlobIndices, GlobIndices_EntType& volumeGlobIndices,
    N_Matrix_Type& diffN_Matrix_nodes,
    N_Matrix_EntType& diffN_Matrix_edges,
    N_Matrix_EntType& diffN_Matrix_faces,
    N_Matrix_EntType& diffN_Matrix_elem) {
  PetscFunctionBegin;
  unsigned int g_dim,nb_Ns;
  switch (fe_ent_ptr->get_ent_type()) {
    case MBTET:
      g_dim = gNTET.size()/4;
      nb_Ns = 4;
      break;
    default:
      SETERRQ(PETSC_COMM_SELF,1,"not implemented yet");
  }
  // nodes
  GlobIndices_Type::iterator nit = nodesGlobIndices.begin();
  for(;nit!=nodesGlobIndices.end();nit++) {
    const MoFEMField* field_ptr = nit->first;
    int rank = field_ptr->get_max_rank();
    int dim = 0,nb_rows = 0;
    switch(field_ptr->get_space()) {
      case H1:
	dim = 3;
	nb_rows = rank*dim;
	break;
      default:
	continue;
    }
    vector<const double*> diff_shape_by_gauss_pt;
    ierr = get_ShapeFunction(NULL,&diff_shape_by_gauss_pt,field_ptr,MBVERTEX); CHKERRQ(ierr);
    assert(diff_shape_by_gauss_pt.size()==g_dim);
    vector< ublas::matrix<FieldData> > &data = diffN_Matrix_nodes[field_ptr];
    data.resize(g_dim);
    if(rank*nb_Ns!=nit->second.size()) SETERRQ(PETSC_COMM_SELF,1,"data inconsitency");
    unsigned int gg = 0;
    for(;gg<g_dim;gg++) {
      ublas::matrix<FieldData> &mat = data[gg];
      mat.resize(nb_rows,rank*nb_Ns);
      mat = ublas::zero_matrix<FieldData>(nb_rows,rank*nb_Ns);
      int rr = 0;
      for(;rr<rank;rr++) {
	int dd = 0;
	for(;dd<dim;dd++) {
	  if(diff_shape_by_gauss_pt[gg] == NULL) SETERRQ(PETSC_COMM_SELF,1,"data inconsitency");
	  ublas::matrix_row<ublas::matrix<FieldData> > mr(mat,rr*dim+dd);
	  for(int nn = 0;nn<4;nn++) {
	    mr(nn*rank + rr) = (diff_shape_by_gauss_pt[gg])[dim*nn+dd];
	  }
	}
      }
      //cerr << mat << endl;
    }
  }
  // edges // faces // volumes
  GlobIndices_EntType* F[] = { &edgesGlobIndices, &facesGlobIndices, &volumeGlobIndices };
  N_Matrix_EntType* FF[] = {  &diffN_Matrix_edges, &diffN_Matrix_faces, &diffN_Matrix_elem };
  for(int ss = 0;ss<3;ss++) {
    for(GlobIndices_EntType::iterator dit = F[ss]->begin();dit!=F[ss]->end();dit++) {
      const MoFEMEntity* ent_ptr = dit->first;
      const MoFEMField* field_ptr = ent_ptr->get_MoFEMField_ptr();
      int rank = field_ptr->get_max_rank();
      int order = ent_ptr->get_max_order();
      int dim = 0,nb_rows = 0;
      switch(field_ptr->get_space()) {
	case H1:
	  dim = 3;
	  nb_rows = rank*dim;
	  break;
	default:
	  continue;
      }
      vector<const double*> diff_shape_by_gauss_pt;
      if(ss<=1) {
	try {
	  int side_number = fe_ent_ptr->get_side_number_ptr(moab,ent_ptr->get_ent())->side_number;
	  ierr = get_ShapeFunction(NULL,&diff_shape_by_gauss_pt,field_ptr,ent_ptr->get_ent_type(),side_number); CHKERRQ(ierr);
	} catch (const char* msg) {
	  SETERRQ(PETSC_COMM_SELF,1,msg);
	}
      } else {
	switch (ent_ptr->get_ent_type()) {
	  case MBTET:
	  case MBPRISM:
	    break;
	  default:
	    SETERRQ(PETSC_COMM_SELF,1,"data inconsitency");
	}
	ierr = get_ShapeFunction(NULL,&diff_shape_by_gauss_pt,field_ptr,ent_ptr->get_ent_type()); CHKERRQ(ierr);
      }
      vector<ublas::matrix<FieldData> > &data = (*FF[ss])[ent_ptr];
      data.resize(g_dim);
      unsigned int nb_dofs = rank*ent_ptr->forder(order); 
      if(nb_dofs!=dit->second.size()) SETERRQ(PETSC_COMM_SELF,1,"data inconsitency");
      unsigned int gg = 0;
      for(;gg<g_dim;gg++) {
	if(diff_shape_by_gauss_pt[gg] == NULL) SETERRQ(PETSC_COMM_SELF,1,"data inconsitency");
	ublas::matrix<FieldData> &mat = data[gg];
	mat.resize(nb_rows,nb_dofs);
	mat = ublas::zero_matrix<FieldData>(nb_rows,nb_dofs);
	for(int rr = 0;rr<rank;rr++) {
	  for(int dd = 0;dd<dim;dd++) {
	    ublas::matrix_row<ublas::matrix<FieldData> > mr(mat,rr*dim+dd);
	    for(int jj = 0;jj<ent_ptr->forder(order);jj++) {
	      mr(jj*rank + rr) =  (diff_shape_by_gauss_pt[gg])[dim*jj + dd];     
	    }
	  }
	}
	//cerr << rank << " " << dim << " " << ent_ptr->forder(order) << endl;
	//cerr << mat << endl;
      }
    }
  }
  PetscFunctionReturn(0);
}
PetscErrorCode FEMethod_LowLevelStudent::GetRowDiffNMatrix_at_GaussPoint() {
  PetscFunctionBegin;
  ierr = GetDiffNMatrix_at_GaussPoint(
    row_nodesGlobIndices,row_edgesGlobIndices,
    row_facesGlobIndices,row_elemGlobIndices,
    row_diffN_Matrix_nodes,row_diffN_Matrix_edges,
    row_diffN_Matrix_faces,row_diffN_Matrix_elem); CHKERRQ(ierr);
  PetscFunctionReturn(0);
}
PetscErrorCode FEMethod_LowLevelStudent::GetColDiffNMatrix_at_GaussPoint() {
  PetscFunctionBegin;
  ierr = GetDiffNMatrix_at_GaussPoint(
    col_nodesGlobIndices,col_edgesGlobIndices,
    col_facesGlobIndices,col_elemGlobIndices,
    col_diffN_Matrix_nodes,col_diffN_Matrix_edges,
    col_diffN_Matrix_faces,col_diffN_Matrix_elem); CHKERRQ(ierr);
  PetscFunctionReturn(0);
}
//TRI
PetscErrorCode FEMethod_LowLevelStudent::ShapeFunctions_TRI(EntityHandle ent,vector<double> &_gNTRI_) {  
  PetscFunctionBegin;
  typedef SideNumber_multiIndex::nth_index<1>::type SideNumber_multiIndex_by_CompositeTag;
  SideNumber_multiIndex_by_CompositeTag& side_table = const_cast<SideNumber_multiIndex_by_CompositeTag&>(fe_ent_ptr->get_side_number_table().get<1>());
  SideNumber* side = fe_ent_ptr->get_side_number_ptr(moab,ent);
  if(side->get_ent_type()!=MBTRI) SETERRQ(PETSC_COMM_SELF,1,"data inconsitency");
  gNTRI = _gNTRI_;
  gNTRI_dim = gNTRI.size()/3;
  if(isH1) {
    switch(fe_ent_ptr->get_ent_type()) {
      case MBTET: {
	//sense and edge order
	const int _faces_edges_[4][3] = { {0,4,3}, {1,5,4}, {3,5,2}, { 0,1,2 } };
	EntityHandle edges[3];
	int _face_edge_sense_[3],_face_edge_offse_[3],_face_edge_side_number_[3];
	int _elem_face_edge_side_number_[3];
	int _face_edge_order_H1[3];
	for(int ee = 0;ee<3;ee++) {
	  _elem_face_edge_side_number_[ee] = _faces_edges_[side->side_number][ee];
	  SideNumber_multiIndex_by_CompositeTag::iterator siit = side_table.find(boost::make_tuple(MBEDGE, _elem_face_edge_side_number_[ee] ));
	  edges[ee] = siit->ent;
	  rval = moab.side_number(side->ent,siit->ent,_face_edge_side_number_[ee],_face_edge_sense_[ee],_face_edge_offse_[ee]); CHKERR_PETSC(rval);
	  assert(_face_edge_side_number_[ee] >= 0); 
	  assert(_face_edge_side_number_[ee] <= 2); 
	  _face_edge_order_H1[ee] = maxOrderEdgeH1[_elem_face_edge_side_number_[ee]];
	} 
	//topology
	const int cannonical_face_sense_p1[4][3] = { {0,1,3}, {1,2,3}, {0,3,2}/**/, {0,2,1}/**/ }; //secon index is offset (positive sense)
	const int cannonical_face_sense_m1[4][3] = { {0,3,1}, {1,3,2}, {0,2,3}, {0,1,2} }; //second index is offset (negative sense)
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
	if(debug) {
	  const EntityHandle *conn_face;
	  int num_nodes;
	  rval = moab.get_connectivity(side->ent,conn_face,num_nodes,true); CHKERR_PETSC(rval);
	  assert(num_nodes==3);
	  assert(conn_face[0] == conn[face_conn[0]]);
	  assert(conn_face[1] == conn[face_conn[1]]);
	  assert(conn_face[2] == conn[face_conn[2]]);
	}
	//nodes
	gNTRIonElem.resize(4);
	assert(side->side_number>=0);
	assert(side->side_number<=3);
	gNTRIonElem[side->side_number].resize(4*gNTRI_dim);
	for(int gg = 0;gg<gNTRI_dim;gg++) {
	  for(int nn = 0;nn<3;nn++) {
	    gNTRIonElem[side->side_number][4*gg + face_conn[nn]] = gNTRI[3*gg + nn];
	  }
	}
	//edges
	H1edgeN_TRI[ent];
	diffH1edgeN_TRI[ent];
	map<EntityHandle,vector<double> >& H1edgeN_TRI_face = H1edgeN_TRI[ent];
	map<EntityHandle,vector<double> >& diffH1edgeN_TRI_face = diffH1edgeN_TRI[ent];
	for(int ee = 0;ee<3;ee++) {
	  H1edgeN_TRI_face[edges[ee]].resize(NBEDGE_H1(_face_edge_order_H1[ee])*gNTRI_dim);
	  diffH1edgeN_TRI_face[edges[ee]].resize(2*NBEDGE_H1(_face_edge_order_H1[ee])*gNTRI_dim);
	}
	double *_edgeN_[] = { &((H1edgeN_TRI_face[edges[0]])[0]), &((H1edgeN_TRI_face[edges[1]])[0]), &((H1edgeN_TRI_face[edges[2]])[0]) };
	double *_diff_edgeN_[] = { &((diffH1edgeN_TRI_face[edges[0]])[0]), &((diffH1edgeN_TRI_face[edges[1]])[0]), &((diffH1edgeN_TRI_face[edges[2]])[0]) };
	ierr = H1_EdgeShapeFunctions_MBTRI(_face_edge_sense_,_face_edge_order_H1,&gNTRI[0],diffNTRI,_edgeN_,_diff_edgeN_,gNTRI_dim); CHKERRQ(ierr);
	//face
	int _face_order_ = maxOrderFaceH1[side->side_number];
	H1faceN_TRI[ent].resize(NBFACE_H1(_face_order_)*gNTRI_dim);
	diffH1faceN_TRI[ent].resize(2*NBFACE_H1(_face_order_)*gNTRI_dim);
	double *_faceN_ = &(H1faceN_TRI[ent][0]);
	double *_diff_faceN_ = &(diffH1faceN_TRI[ent][0]);
	ierr = H1_FaceShapeFunctions_MBTRI(_face_order_,&gNTRI[0],diffNTRI,_faceN_,_diff_faceN_,gNTRI_dim); CHKERRQ(ierr);
      }
      break;
      default:
	SETERRQ(PETSC_COMM_SELF,1,"not implemented");
    }
  }
  PetscFunctionReturn(0);
}
PetscErrorCode FEMethod_LowLevelStudent::GetNMatrix_at_FaceGaussPoint(
    EntityHandle ent,const string& field_name,
    GlobIndices_Type& nodesGlobIndices, 
    GlobIndices_EntType& edgesGlobIndices,
    GlobIndices_EntType& facesGlobIndices,
    N_Matrix_Type& N_Matrix_nodes,
    N_Matrix_EntType& N_Matrix_edges,
    N_Matrix_EntType& N_Matrix_faces,
    EntityType type,EntityHandle edge_handle) {
  PetscFunctionBegin;
  SideNumber* side = fe_ent_ptr->get_side_number_ptr(moab,ent);
  if(side->get_ent_type()!=MBTRI) SETERRQ(PETSC_COMM_SELF,1,"entity has to be face of type MBTRI");
  unsigned int g_dim,nb_Ns;
  switch (fe_ent_ptr->get_ent_type()) {
    case MBTET:
      g_dim = get_dim_gNTRI();
      nb_Ns = 4;
      break;
    default:
      SETERRQ(PETSC_COMM_SELF,1,"not implemented yet");
  }
  int side_number = side->side_number;
  // nodes
  if((type == MBVERTEX)||(type == MBMAXTYPE)) {
    GlobIndices_Type::iterator nit = nodesGlobIndices.begin();
    for(;nit!=nodesGlobIndices.end();nit++) {
      const MoFEMField* field_ptr = nit->first;
      if(field_ptr->get_name()!=field_name) continue;
      vector<double>& gNTRIonELEM = gNTRIonElem[side_number];
      if(nb_Ns*g_dim != gNTRIonELEM.size()) SETERRQ(PETSC_COMM_SELF,1,"data inconsitency");
      vector< ublas::matrix<FieldData> > &data = N_Matrix_nodes[field_ptr];
      data.resize(g_dim);
      int rank = field_ptr->get_max_rank();
      if(rank*nb_Ns!=nit->second.size()) SETERRQ(PETSC_COMM_SELF,1,"data inconsitency");
      unsigned int gg = 0;
      for(;gg<g_dim;gg++) {
	ublas::matrix<FieldData> &mat = data[gg];
	mat.resize(rank,rank*nb_Ns);
	mat = ublas::zero_matrix<FieldData>(rank,rank*nb_Ns);
	int rr = 0;
	for(;rr<rank;rr++) {
	  ublas::matrix_row<ublas::matrix<FieldData> > mr(mat,rr);
	  for(unsigned int jj = 0;jj<nb_Ns;jj++) mr(rank*jj + rr) = gNTRIonELEM[gg*nb_Ns + jj];
	}
      }
    }
  }
  //edges
  if((type == MBEDGE)||(type == MBMAXTYPE)) {
    map<EntityHandle,vector<double> > &H1edgeN_TRI_face = H1edgeN_TRI[ent];
    for(GlobIndices_EntType::iterator eiit = edgesGlobIndices.begin();eiit!=edgesGlobIndices.end();eiit++) {
      const MoFEMEntity* ent_ptr = eiit->first;
      const MoFEMField* field_ptr = ent_ptr->get_MoFEMField_ptr();
      if(field_ptr->get_name()!=field_name) continue;
      EntityHandle edge = eiit->first->get_ent();
      if(edge_handle!=no_handle) {
	if(edge_handle!=ent_ptr->get_ent()) continue;
      }
      map<EntityHandle,vector<double> >::iterator mit = H1edgeN_TRI_face.find(edge);
      vector<ublas::matrix<FieldData> > &data = N_Matrix_edges[ent_ptr];
      data.resize(g_dim);
      int rank = field_ptr->get_max_rank();
      int order = ent_ptr->get_max_order();
      unsigned int nb_dofs = rank*ent_ptr->forder(order); 
      if(mit!=H1edgeN_TRI_face.end()) {
	if(ent_ptr->forder(order)*g_dim != mit->second.size()) 
	  SETERRQ(PETSC_COMM_SELF,1,"data inconsitency");
      }
      if(nb_dofs!=eiit->second.size()) SETERRQ(PETSC_COMM_SELF,1,"data inconsitency");
      unsigned int gg = 0;
      for(;gg<g_dim;gg++) {
	ublas::matrix<FieldData> &mat = data[gg];
	mat.resize(rank,nb_dofs);
	mat = ublas::zero_matrix<FieldData>(rank,nb_dofs);
	if(mit==H1edgeN_TRI_face.end()) continue;
	for(int rr = 0;rr<rank;rr++) {
	  ublas::matrix_row<ublas::matrix<FieldData> > mr(mat,rr);
	  for(int jj = 0;jj<ent_ptr->forder(order);jj++) {
	    mr(rank*jj + rr) = (mit->second)[gg*ent_ptr->forder(order)+jj];
	  }
	}
      }
    }
  }
  //faces
  if((type == MBTRI)||(type == MBMAXTYPE)) {
    vector<double> &H1faceN_TRI_face = H1faceN_TRI[ent];
    for(GlobIndices_EntType::iterator fiit = facesGlobIndices.begin();fiit!=facesGlobIndices.end();fiit++) {
      const MoFEMEntity* ent_ptr = fiit->first;
      const MoFEMField* field_ptr = ent_ptr->get_MoFEMField_ptr();
      if(field_ptr->get_name()!=field_name) continue;
      EntityHandle face = fiit->first->get_ent();
      vector<ublas::matrix<FieldData> > &data = N_Matrix_faces[ent_ptr];
      data.resize(g_dim);
      int rank = field_ptr->get_max_rank();
      int order = ent_ptr->get_max_order();
      unsigned int nb_dofs = rank*ent_ptr->forder(order); 
      if(ent==face) {
	if(ent_ptr->forder(order)*g_dim != H1faceN_TRI_face.size()) 
	  SETERRQ1(PETSC_COMM_SELF,1,"data inconsitency (side_number = %u)",side->side_number);
      }
      if(nb_dofs*g_dim!=rank*H1faceN_TRI_face.size()) SETERRQ1(PETSC_COMM_SELF,1,"data inconsitency (side_number = %u)",side->side_number);
      unsigned int gg = 0;
      for(;gg<g_dim;gg++) {
	ublas::matrix<FieldData> &mat = data[gg];
	mat.resize(rank,nb_dofs);
	mat = ublas::zero_matrix<FieldData>(rank,nb_dofs);
	if(ent!=face) continue;
	for(int rr = 0;rr<rank;rr++) {
	  ublas::matrix_row<ublas::matrix<FieldData> > mr(mat,rr);
	  for(int jj = 0;jj<ent_ptr->forder(order);jj++) {
	    mr(rank*jj + rr) = H1faceN_TRI_face[gg*ent_ptr->forder(order)+jj];
	  }
	}
      }
    }
  }
  PetscFunctionReturn(0);
}


}
