/* Copyright (C) 2013, Lukasz Kaczmarczyk (likask AT wp.pl)
 * --------------------------------------------------------------
 * FIXME: DESCRIPTION
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

#include "FEMethod_ComplexForLazy.hpp"
#include "FEM.h"
#include "complex_for_lazy.h"

extern "C" {
void tetcircumcenter_tp(double a[3],double b[3],double c[3], double d[3],
  double circumcenter[3],double *xi,double *eta,double *zeta);
void tricircumcenter3d_tp(double a[3],double b[3],double c[3],
  double circumcenter[3],double *xi,double *eta);
}

namespace MoFEM {

FEMethod_ComplexForLazy::FEMethod_ComplexForLazy(FieldInterface& _mField,BaseDirihletBC *_dirihlet_bc_method_ptr,
    analysis _type,
    double _lambda,double _mu, int _verbose): 
    FEMethod_ComplexForLazy_Data(_mField,_dirihlet_bc_method_ptr,_verbose), 
    type_of_analysis(_type),type_of_forces(conservative),lambda(_lambda),mu(_mu),eps(1e-10),
    spatial_field_name("SPATIAL_POSITION"),
    material_field_name("MESH_NODE_POSITIONS") {
  order_edges.resize(6);
  order_faces.resize(4);
  edgeNinvJac.resize(6);
  faceNinvJac.resize(4);
  diff_edgeNinvJac.resize(6);
  diff_faceNinvJac.resize(4);
  //Tangent_HH_hierachical
  KedgeH_data.resize(6);
  KfaceH_data.resize(4);
  //Tangent_hh_hierachical
  Kedgeh_data.resize(6);
  Kfaceh_data.resize(4);
  //Tangent_hh_hierachical_edge
  Khedge_data.resize(6);
  KHedge_data.resize(6);
  Khh_volumeedge_data.resize(6);
  Khh_edgeedge_data.resize(6,6);
  Khh_faceedge_data.resize(4,6);
  //Tangent_hh_hierachical_face
  Khface_data.resize(6);
  KHface_data.resize(6);
  Khh_volumeface_data.resize(4);
  Khh_faceface_data.resize(4,4);
  Khh_edgeface_data.resize(6,4);
  //Tangent_hh_hierachical_volume
  Khh_edgevolume_data.resize(6);
  Khh_facevolume_data.resize(4);
  //
  Fint_h.resize(12);
  Fint_h_edge_data.resize(6);
  Fint_h_face_data.resize(4);
  //
  Fint_H.resize(12);
  //
  g_NTET.resize(4*45);
  ShapeMBTET(&g_NTET[0],G_TET_X45,G_TET_Y45,G_TET_Z45,45);
  g_TET_W = G_TET_W45;
  //
  g_NTRI.resize(3*13);
  ShapeMBTRI(&g_NTRI[0],G_TRI_X13,G_TRI_Y13,13); 
  g_TRI_dim = 13;
  g_TRI_W = G_TRI_W13;

  propeties_from_BlockSet_Mat_ElasticSet = false;
  for(_IT_CUBITMESHSETS_BY_BCDATA_TYPE_FOR_LOOP_(mField,BlockSet|Mat_ElasticSet,it)) {
    propeties_from_BlockSet_Mat_ElasticSet = true;
  }

}
PetscErrorCode FEMethod_ComplexForLazy::GetMatParameters(double *_lambda,double *_mu,void *ptr_matctx) {
  PetscFunctionBegin;

  *_lambda = lambda;
  *_mu = mu;

  if(propeties_from_BlockSet_Mat_ElasticSet) {
    EntityHandle ent = fe_ptr->get_ent();
    for(_IT_CUBITMESHSETS_BY_BCDATA_TYPE_FOR_LOOP_(mField,BlockSet|Mat_ElasticSet,it)) {

      Mat_Elastic mydata;
      ierr = it->get_attribute_data_structure(mydata); CHKERRQ(ierr);

      Range meshsets;
      rval = moab.get_entities_by_type(it->meshset,MBENTITYSET,meshsets,true); CHKERR_PETSC(rval);
      meshsets.insert(it->meshset);
      for(Range::iterator mit = meshsets.begin();mit != meshsets.end(); mit++) {
	if( moab.contains_entities(*mit,&ent,1) ) {
	  *_lambda = LAMBDA(mydata.data.Young,mydata.data.Poisson);
	  *_mu = MU(mydata.data.Young,mydata.data.Poisson);
	  PetscFunctionReturn(0);  
	}
      }

    }

    SETERRQ(PETSC_COMM_SELF,1,
      "Element is not in elastic block, however you run non-linear elastic analysis with that element\n"
      "top tip: check if you update block sets after mesh refinments or interface insertion");

  }

  PetscFunctionReturn(0);
}
PetscErrorCode FEMethod_ComplexForLazy::OpComplexForLazyStart() {
  PetscFunctionBegin;
  ierr = OpStudentStart_TET(g_NTET); CHKERRQ(ierr);
  if(mesh_quality_analysis&type_of_analysis) {
    ierr = PetscOptionsGetReal("","-my_alpha2",&alpha2,&flg_alpha2); CHKERRQ(ierr);
    ierr = PetscOptionsGetReal("","-my_gamma",&gamma,&flg_gamma); CHKERRQ(ierr);
    double def_quality = -1;
    rval = moab.tag_get_handle("QUALITY0",1,MB_TYPE_DOUBLE,th_quality0,MB_TAG_CREAT|MB_TAG_SPARSE,&def_quality); 
    if(rval==MB_ALREADY_ALLOCATED) rval = MB_SUCCESS;
    CHKERR(rval);
    rval = moab.tag_get_handle("QUALITY",1,MB_TYPE_DOUBLE,th_quality,MB_TAG_CREAT|MB_TAG_SPARSE,&def_quality); 
    if(rval==MB_ALREADY_ALLOCATED) rval = MB_SUCCESS;
    CHKERR(rval);
    rval = moab.tag_get_handle("B",1,MB_TYPE_DOUBLE,th_b,MB_TAG_CREAT|MB_TAG_SPARSE,&def_quality); 
    if(rval==MB_ALREADY_ALLOCATED) rval = MB_SUCCESS;
    CHKERR(rval);
    EntityHandle ent = fe_ptr->get_ent();
    rval = moab.tag_get_by_ptr(th_quality0,&ent,1,(const void **)&quality0); CHKERR(rval);
    rval = moab.tag_get_by_ptr(th_quality,&ent,1,(const void **)&quality); CHKERR(rval);
    rval = moab.tag_get_by_ptr(th_b,&ent,1,(const void **)&b); CHKERR(rval);
  }
  PetscFunctionReturn(0);
}
PetscErrorCode FEMethod_ComplexForLazy::GetIndicesRow(
  vector<vector<DofIdx> > &RowGlob,string &field_name) {
  PetscFunctionBegin;
  try{
  switch (fe_ent_ptr->get_ent_type()) {
    case MBTET: {
      RowGlob.resize(1+6+4+1);
      try {
      //nodes
      ierr = GetRowGlobalIndices(field_name,RowGlob[0]); CHKERRQ(ierr);
      //edges
      int ee = 0;
      for(;ee<6;ee++) { //edges matrices
	ierr = GetRowGlobalIndices(field_name,MBEDGE,RowGlob[1+ee],ee); CHKERRQ(ierr);
      }
      assert(ee == 6);
      //faces
      int ff = 0;
      for(;ff<4;ff++) { //faces matrices
	ierr = GetRowGlobalIndices(field_name,MBTRI,RowGlob[1+ee+ff],ff); CHKERRQ(ierr);
      }
      assert(ff == 4);
      //volumes
      ierr = GetRowGlobalIndices(field_name,MBTET,RowGlob[1+ee+ff]); CHKERRQ(ierr);
      } catch (const char* msg) {
	SETERRQ(PETSC_COMM_SELF,1,msg);
      } 
    }
    break;
    default:
      SETERRQ(PETSC_COMM_SELF,1,"sorry.. I don't know what to do");
  }
  } catch (const std::exception& ex) {
      ostringstream ss;
      ss << "throw in method: " << ex.what() << endl;
      SETERRQ(PETSC_COMM_SELF,1,ss.str().c_str());
  } 
  PetscFunctionReturn(0);
}
PetscErrorCode FEMethod_ComplexForLazy::GetIndicesCol(
  vector<vector<DofIdx> > &ColGlob,string &field_name) {
  PetscFunctionBegin;
  try{
  switch (fe_ent_ptr->get_ent_type()) {
    case MBTET: {
      ColGlob.resize(1+6+4+1);
      try {
      //nodes
      ierr = GetColGlobalIndices(field_name,ColGlob[0]); CHKERRQ(ierr);
      //edges
      int ee = 0;
      for(;ee<6;ee++) { //edges matrices
	ierr = GetColGlobalIndices(field_name,MBEDGE,ColGlob[1+ee],ee); CHKERRQ(ierr);
      }
      assert(ee == 6);
      //faces
      int ff = 0;
      for(;ff<4;ff++) { //faces matrices
	ierr = GetColGlobalIndices(field_name,MBTRI,ColGlob[1+ee+ff],ff); CHKERRQ(ierr);
      }
      assert(ff == 4);
      //volumes
      ierr = GetColGlobalIndices(field_name,MBTET,ColGlob[1+ee+ff]); CHKERRQ(ierr);
      } catch (const char* msg) {
	SETERRQ(PETSC_COMM_SELF,1,msg);
      } 
    }
    break;
    default:
      SETERRQ(PETSC_COMM_SELF,1,"sorry.. I don't know what to do");
  }
  } catch (const std::exception& ex) {
      ostringstream ss;
      ss << "throw in method: " << ex.what() << endl;
      SETERRQ(PETSC_COMM_SELF,1,ss.str().c_str());
  } 
  PetscFunctionReturn(0);
}
PetscErrorCode FEMethod_ComplexForLazy::GetIndices(
  vector<vector<DofIdx> > &RowGlob,vector<vector<DofIdx> > &ColGlob,
  string &field_name) {
  PetscFunctionBegin;
  ierr = GetIndicesRow(RowGlob,field_name); CHKERRQ(ierr);
  ierr = GetIndicesCol(ColGlob,field_name); CHKERRQ(ierr);
  PetscFunctionReturn(0);
}
PetscErrorCode FEMethod_ComplexForLazy::GetData(
  vector<ublas::vector<double> >& dofs_edge_data,vector<double*>& dofs_edge,
  vector<ublas::vector<double> >& dofs_face_data,vector<double*>& dofs_face,
  ublas::vector<double>& dofs_volume,ublas::vector<double>& dofs_nodes,
  string &field_name) {
  PetscFunctionBegin;
  try{
  switch (fe_ent_ptr->get_ent_type()) {
    case MBTET: {
      dofs_nodes.resize(12);
      dofs_edge_data.resize(6);
      dofs_face_data.resize(4);
      dofs_edge.resize(6);
      dofs_face.resize(4);
      try {
      //data edge
      int ee = 0;
      for(;ee<6;ee++) {
	FEDofMoFEMEntity_multiIndex::index<Composite_Name_Type_And_Side_Number_mi_tag>::type::iterator eiit,hi_eiit;
	eiit = data_multiIndex->get<Composite_Name_Type_And_Side_Number_mi_tag>().lower_bound(boost::make_tuple(field_name,MBEDGE,ee));
	hi_eiit = data_multiIndex->get<Composite_Name_Type_And_Side_Number_mi_tag>().upper_bound(boost::make_tuple(field_name,MBEDGE,ee));
	if(eiit!=hi_eiit) {
	  assert(eiit->side_number_ptr->side_number==ee);
	  dofs_edge_data[ee].resize(distance(eiit,hi_eiit));
	  order_edges[ee] = eiit->get_max_order();
	  dofs_edge[ee] = &dofs_edge_data[ee].data()[0];
	  assert(dofs_edge_data[ee].size() == 3*(unsigned int)NBEDGE_H1(order_edges[ee]));
	  map<UId,FieldData> map_edge_dofs;
	  for(;eiit!=hi_eiit;eiit++) dofs_edge_data[ee][eiit->get_EntDofIdx()] = eiit->get_FieldData(); 
	} else {
	  order_edges[ee] = 0;
	  dofs_edge[ee] = NULL;
	}
      }
      //data face
      int ff = 0;
      for(;ff<4;ff++) {
	FEDofMoFEMEntity_multiIndex::index<Composite_Name_Type_And_Side_Number_mi_tag>::type::iterator fiit,hi_fiit;
	fiit = data_multiIndex->get<Composite_Name_Type_And_Side_Number_mi_tag>().lower_bound(boost::make_tuple(field_name,MBTRI,ff));
	hi_fiit = data_multiIndex->get<Composite_Name_Type_And_Side_Number_mi_tag>().upper_bound(boost::make_tuple(field_name,MBTRI,ff));
	if(fiit!=hi_fiit) {
	  dofs_face_data[ff].resize(distance(fiit,hi_fiit));
	  order_faces[ff] = fiit->get_max_order();
	  dofs_face[ff] = &dofs_face_data[ff].data()[0];
	  assert(dofs_face_data[ff].size() == 3*(unsigned int)NBFACE_H1(order_faces[ff]));
	  for(;fiit!=hi_fiit;fiit++) dofs_face_data[ff][fiit->get_EntDofIdx()] = fiit->get_FieldData(); 
	} else {
	  order_faces[ff] = 0;
	  dofs_face[ff] = 0;
	}
      }
      //data voolume
      FEDofMoFEMEntity_multiIndex::index<Composite_mi_tag2>::type::iterator viit,hi_viit;
      viit = data_multiIndex->get<Composite_mi_tag2>().lower_bound(boost::make_tuple(field_name,MBTET));
      hi_viit = data_multiIndex->get<Composite_mi_tag2>().upper_bound(boost::make_tuple(field_name,MBTET));
      if(viit!=hi_viit) {
	order_volume = viit->get_max_order();
	dofs_volume.resize(distance(viit,hi_viit));
	assert(dofs_volume.size() == (unsigned int)3*NBVOLUME_H1(order_volume));
	for(;viit!=hi_viit;viit++) dofs_volume[viit->get_EntDofIdx()] = viit->get_FieldData(); 
      } else {
	order_volume = 0;
      }
      //data nodes
      FEDofMoFEMEntity_multiIndex::index<Composite_Name_Type_And_Side_Number_mi_tag>::type::iterator niit,hi_niit;
      niit = data_multiIndex->get<Composite_Name_Type_And_Side_Number_mi_tag>().lower_bound(boost::make_tuple(field_name,MBVERTEX,0));
      hi_niit = data_multiIndex->get<Composite_Name_Type_And_Side_Number_mi_tag>().upper_bound(boost::make_tuple(field_name,MBVERTEX,4));
      if(distance(niit,hi_niit)!=12) SETERRQ(PETSC_COMM_SELF,1,"I can not find dofs on vertices, it should be 12 dofs (i.e. 4 nodes and 3 dofs for each node)");
      for(int dd = 0;niit!=hi_niit;niit++,dd++) {
	dofs_nodes[3*niit->side_number_ptr->side_number+niit->get_EntDofIdx()] = niit->get_FieldData(); 
      }
      } catch (const char* msg) {
	SETERRQ(PETSC_COMM_SELF,1,msg);
      } 
    }  break;
    default:
      SETERRQ(PETSC_COMM_SELF,1,"sorry.. I don't know what to do");
  }
  } catch (const std::exception& ex) {
      ostringstream ss;
      ss << "throw in method: " << ex.what() << endl;
      SETERRQ(PETSC_COMM_SELF,1,ss.str().c_str());
  } 
  PetscFunctionReturn(0);
}

PetscErrorCode FEMethod_ComplexForLazy::GetIndicesSpatial() {
  PetscFunctionBegin;
  ierr = GetIndices(RowGlobSpatial,ColGlobSpatial,spatial_field_name); CHKERRQ(ierr);
  ierr = GetData(dofs_x_edge_data,dofs_x_edge,
    dofs_x_face_data,dofs_x_face,
    dofs_x_volume,dofs_x,
    spatial_field_name); CHKERRQ(ierr);
  PetscFunctionReturn(0);
}
PetscErrorCode FEMethod_ComplexForLazy::GetIndicesMaterial() {
  PetscFunctionBegin;
  ierr = GetIndices(RowGlobMaterial,ColGlobMaterial,material_field_name); CHKERRQ(ierr);
  PetscFunctionReturn(0);
}
PetscErrorCode FEMethod_ComplexForLazy::GetDofs_X_FromElementData() {
  PetscFunctionBegin;
  dofs_X.resize(12);
  copy(coords.begin(),coords.end(),dofs_X.begin());
  FEDofMoFEMEntity_multiIndex::index<Composite_Name_Type_And_Side_Number_mi_tag>::type::iterator niit,hi_niit;
  niit = data_multiIndex->get<Composite_Name_Type_And_Side_Number_mi_tag>().lower_bound(boost::make_tuple(material_field_name,MBVERTEX,0));
  hi_niit = data_multiIndex->get<Composite_Name_Type_And_Side_Number_mi_tag>().upper_bound(boost::make_tuple(material_field_name,MBVERTEX,4));
  for(;niit!=hi_niit;niit++) {
    dofs_X[3*niit->side_number_ptr->side_number+niit->get_EntDofIdx()] = niit->get_FieldData();
  }
  PetscFunctionReturn(0);
}
PetscErrorCode FEMethod_ComplexForLazy::GetTangent() {
  PetscFunctionBegin;
  try {
  switch (fe_ent_ptr->get_ent_type()) {
  case MBTET: {
    int ee = 0;
    for(;ee<6;ee++) {
      diff_edgeNinvJac[ee] = &*diffH1edgeNinvJac[ee].begin(); 
    }
    int ff = 0;
    for(;ff<4;ff++) {
      diff_faceNinvJac[ff] = &*diffH1faceNinvJac[ff].begin();
    }
    diff_volumeNinvJac = &*diffH1elemNinvJac.begin();
    double center[3]; 
    tetcircumcenter_tp(&coords[0],&coords[3],&coords[6],&coords[9],center,NULL,NULL,NULL); 
    cblas_daxpy(3,-1,&coords[0],1,center,1);
    double r = cblas_dnrm2(3,center,1);
    int g_dim = get_dim_gNTET();
    if(type_of_analysis&spatail_analysis) {
	assert(12 == RowGlobSpatial[0].size());
	KHh.resize(12,12);
	Khh.resize(12,12);
	ee = 0;
	for(;ee<6;ee++) {
	  assert(3*(unsigned int)NBEDGE_H1(order_edges[ee]) == RowGlobSpatial[1+ee].size());
	  Kedgeh_data[ee].resize(RowGlobSpatial[1+ee].size(),12);
	  Kedgeh[ee] = &*Kedgeh_data[ee].data().begin();
	  KHedge_data[ee].resize(12,ColGlobSpatial[1+ee].size());
	  KHedge[ee] = &*KHedge_data[ee].data().begin();
	  Khedge_data[ee].resize(12,ColGlobSpatial[1+ee].size());
	  Khedge[ee] = &*Khedge_data[ee].data().begin();
	  for(int eee = 0;eee<6;eee++) {
	    Khh_edgeedge_data(eee,ee).resize(RowGlobSpatial[1+eee].size(),ColGlobSpatial[1+ee].size());
	    Khh_edgeedge[eee][ee] = &*Khh_edgeedge_data(eee,ee).data().begin();
	  }
	  for(int fff = 0;fff<4;fff++) {
	    Khh_faceedge_data(fff,ee).resize(RowGlobSpatial[1+6+fff].size(),ColGlobSpatial[1+ee].size());
	    Khh_faceedge[fff][ee] = &*Khh_faceedge_data(fff,ee).data().begin();
	  }
	  Khh_volumeedge_data[ee].resize(RowGlobSpatial[i_volume].size(),ColGlobSpatial[1+ee].size());
	  Khh_volumeedge[ee] = &*Khh_volumeedge_data[ee].data().begin();
	  Khh_edgevolume_data[ee].resize(RowGlobSpatial[1+ee].size(),ColGlobSpatial[i_volume].size());
	  Khh_edgevolume[ee] = &*Khh_edgevolume_data[ee].data().begin();
	}
	ff = 0;
	for(;ff<4;ff++) {
	  assert(3*(unsigned int)NBFACE_H1(order_faces[ff]) == RowGlobSpatial[1+6+ff].size());
	  Kfaceh_data[ff].resize(RowGlobSpatial[1+6+ff].size(),12);
	  Kfaceh[ff] = &*Kfaceh_data[ff].data().begin();
	  Khface_data[ff].resize(12,ColGlobSpatial[1+6+ff].size());
	  Khface[ff] = &*Khface_data[ff].data().begin();
	  KHface_data[ff].resize(12,ColGlobSpatial[1+6+ff].size());
	  KHface[ff] = &*KHface_data[ff].data().begin();
	  for(int fff = 0;fff<4;fff++) {
	    Khh_faceface_data(fff,ff).resize(RowGlobSpatial[1+6+fff].size(),ColGlobSpatial[1+6+ff].size());
	    Khh_faceface[fff][ff] = &*Khh_faceface_data(fff,ff).data().begin();
	  }
	  for(int eee = 0;eee<6;eee++) {
	    Khh_edgeface_data(eee,ff).resize(RowGlobSpatial[1+eee].size(),ColGlobSpatial[1+6+ff].size());
	    Khh_edgeface[eee][ff] = &*Khh_edgeface_data(eee,ff).data().begin();
	  }
	  Khh_volumeface_data[ff].resize(RowGlobSpatial[i_volume].size(),ColGlobSpatial[1+6+ff].size());
	  Khh_volumeface[ff] = &*Khh_volumeface_data[ff].data().begin();
	  Khh_facevolume_data[ff].resize(RowGlobSpatial[1+6+ff].size(),ColGlobSpatial[i_volume].size());
	  Khh_facevolume[ff] = &*Khh_facevolume_data[ff].data().begin();
	}
	assert(3*(unsigned int)NBVOLUME_H1(order_volume) == RowGlobSpatial[i_volume].size());
	Kvolumeh.resize(RowGlobSpatial[i_volume].size(),12);
	KHvolume.resize(12,ColGlobSpatial[i_volume].size());
	Khvolume.resize(12,ColGlobSpatial[i_volume].size());
	Khh_volumevolume.resize(RowGlobSpatial[i_volume].size(),ColGlobSpatial[i_volume].size());
    }
    if(type_of_analysis&material_analysis) {
      assert(12 == RowGlobMaterial[0].size());
      KHH.resize(12,12);
      KhH.resize(12,12);
      ee = 0;
      for(;ee<6;ee++) {
	FEDofMoFEMEntity_multiIndex::index<Composite_Name_Type_And_Side_Number_mi_tag>::type::iterator eiit;
	eiit = data_multiIndex->get<Composite_Name_Type_And_Side_Number_mi_tag>().find(boost::make_tuple(spatial_field_name,MBEDGE,ee));
	if(eiit==data_multiIndex->get<Composite_Name_Type_And_Side_Number_mi_tag>().end()) {
	  order_edges[ee] = 0;
	} else {
	  order_edges[ee] = eiit->get_max_order();
	}
	if(3*(unsigned int)NBEDGE_H1(order_edges[ee]) != dofs_x_edge_data[ee].size()) {
	  SETERRQ2(PETSC_COMM_SELF,1,"data inconsistency 3*NBEDGE_H1(order_edges[ee]) = %d != dofs_x_edge_data[ee].size() = %d ",
	    3*NBEDGE_H1(order_edges[ee]),dofs_x_edge_data[ee].size());
	}
	KedgeH_data[ee].resize(dofs_x_edge_data[ee].size(),12);
	KedgeH[ee] = &*KedgeH_data[ee].data().begin();
      }
      ff = 0;
      for(;ff<4;ff++) {
	FEDofMoFEMEntity_multiIndex::index<Composite_Name_Type_And_Side_Number_mi_tag>::type::iterator fiit;
	fiit = data_multiIndex->get<Composite_Name_Type_And_Side_Number_mi_tag>().find(boost::make_tuple(spatial_field_name,MBTRI,ff));
	if(fiit==data_multiIndex->get<Composite_Name_Type_And_Side_Number_mi_tag>().end()) {
	  order_faces[ff] = 0; 
	} else {
	  order_faces[ff] = fiit->get_max_order();
	}
	assert(3*(unsigned int)NBFACE_H1(order_faces[ff]) == dofs_x_face_data[ff].size());
	KfaceH_data[ff].resize(dofs_x_face_data[ff].size(),12);
	KfaceH[ff] = &*KfaceH_data[ff].data().begin();
      }
      FEDofMoFEMEntity_multiIndex::index<Composite_mi_tag2>::type::iterator viit;
      viit = data_multiIndex->get<Composite_mi_tag2>().find(boost::make_tuple(spatial_field_name,MBTET));
      if(viit==data_multiIndex->get<Composite_mi_tag2>().end()) {
	order_volume = 0;
      } else {
	order_volume = viit->get_max_order();
      }
      KvolumeH.resize(dofs_x_volume.size(),12);
    }
    ierr = GetDofs_X_FromElementData(); CHKERRQ(ierr);
    unsigned int sub_analysis_type[3] = { spatail_analysis, material_analysis, mesh_quality_analysis };
    double _lambda,_mu;
    if( (spatail_analysis|material_analysis)&type_of_analysis ) {
      ierr = GetMatParameters(&_lambda,&_mu,ptr_matctx); CHKERRQ(ierr);
    }
    for(int ss = 0;ss<3;ss++) {
      switch(sub_analysis_type[ss]&type_of_analysis) {
	case spatail_analysis: {
	  ierr = Tangent_hh_hierachical(&order_edges[0],&order_faces[0],order_volume,V,eps*r,_lambda,_mu,ptr_matctx,
	    &diffNTETinvJac[0],&diff_edgeNinvJac[0],&diff_faceNinvJac[0],&diff_volumeNinvJac[0], 
	    &dofs_X.data()[0],&dofs_x[0],&dofs_x_edge[0],&dofs_x_face[0],&*dofs_x_volume.data().begin(), 
	    &*Khh.data().begin(),&*KHh.data().begin(),Kedgeh,Kfaceh,&*Kvolumeh.data().begin(),
	    g_dim,g_TET_W); CHKERRQ(ierr);
	  ierr = Tangent_hh_hierachical_edge(&order_edges[0],&order_faces[0],order_volume,V,eps*r,_lambda,_mu,ptr_matctx, 
	    &diffNTETinvJac[0],&diff_edgeNinvJac[0],&diff_faceNinvJac[0],&diff_volumeNinvJac[0], 
	    &dofs_X.data()[0],&dofs_x[0],&dofs_x_edge[0],&dofs_x_face[0],&*dofs_x_volume.data().begin(), 
	    &Khedge[0],&KHedge[0],Khh_edgeedge,Khh_faceedge,Khh_volumeedge, 
	    g_dim,g_TET_W); CHKERRQ(ierr);
	  ierr = Tangent_hh_hierachical_face(&order_edges[0],&order_faces[0],order_volume,V,eps*r,_lambda,_mu,ptr_matctx, 
	    &diffNTETinvJac[0],&diff_edgeNinvJac[0],&diff_faceNinvJac[0],&diff_volumeNinvJac[0], 
	    &dofs_X.data()[0],&dofs_x[0],&dofs_x_edge[0],&dofs_x_face[0],&*dofs_x_volume.data().begin(), 
	    &Khface[0],&KHface[0],Khh_edgeface,Khh_faceface,Khh_volumeface, 
	    g_dim,g_TET_W); CHKERRQ(ierr);
	  ierr = Tangent_hh_hierachical_volume(&order_edges[0],&order_faces[0],order_volume,V,eps*r,_lambda,_mu,ptr_matctx, 
	    &diffNTETinvJac[0],&diff_edgeNinvJac[0],&diff_faceNinvJac[0],diff_volumeNinvJac, 
	    &dofs_X.data()[0],&dofs_x[0],&dofs_x_edge[0],&dofs_x_face[0],&*dofs_x_volume.data().begin(), 
	    &*Khvolume.data().begin(),&*KHvolume.data().begin(),Khh_edgevolume,Khh_facevolume,&*Khh_volumevolume.data().begin(), 
	    g_dim,g_TET_W); CHKERRQ(ierr);
	}
	break;
	case material_analysis: {
	ierr = Tangent_HH_hierachical(&order_edges[0],&order_faces[0],order_volume,V,eps*r,_lambda,_mu,ptr_matctx,
	  &diffNTETinvJac[0],&diff_edgeNinvJac[0],&diff_faceNinvJac[0],&diff_volumeNinvJac[0], 
	  &dofs_X.data()[0],&dofs_x[0],&dofs_x_edge[0],&dofs_x_face[0],&*dofs_x_volume.data().begin(), 
	  &*KHH.data().begin(),&*KhH.data().begin(),KedgeH,KfaceH,&*KvolumeH.data().begin(),
	  g_dim,g_TET_W); CHKERRQ(ierr);
	}
	break;
	case mesh_quality_analysis: {
	  if(!flg_alpha2) SETERRQ(PETSC_COMM_SELF,1,"-my_alpha2 is not set");
	  if(!flg_gamma) SETERRQ(PETSC_COMM_SELF,1,"-my_gamma is not set");
	  KHH.resize(12,12);
	  double coords_edges[2*3*6]; 
	  double alpha2_array[4] = { alpha2,alpha2,alpha2,alpha2 }; 
	  ierr = get_edges_from_elem_coords(&coords[0],coords_edges); CHKERRQ(ierr);
	  ierr = quality_volume_length_K(eps*r,V,alpha2_array,gamma,&diffNTETinvJac[0],coords_edges,&dofs_X.data()[0],NULL,&*KHH.data().begin(),NULL); CHKERRQ(ierr);
	}
	default:
	  continue;
      }
    }
  }
  break;
  default:
    SETERRQ(PETSC_COMM_SELF,1,"not implemented");
  }
  } catch (const std::exception& ex) {
    ostringstream ss;
    ss << "throw in method: " << ex.what() << endl;
    SETERRQ(PETSC_COMM_SELF,1,ss.str().c_str());
  } 
  PetscFunctionReturn(0);
}
PetscErrorCode FEMethod_ComplexForLazy::get_edges_from_elem_coords(double *coords,double *coords_edges) {
  PetscFunctionBegin;
  cblas_dcopy(3,&coords[0*3],1,&coords_edges[0* 3*2+0],1); 
  cblas_dcopy(3,&coords[1*3],1,&coords_edges[0* 3*2+3],1); 
  cblas_dcopy(3,&coords[0*3],1,&coords_edges[1* 3*2+0],1); 
  cblas_dcopy(3,&coords[2*3],1,&coords_edges[1* 3*2+3],1); 
  cblas_dcopy(3,&coords[0*3],1,&coords_edges[2* 3*2+0],1);
  cblas_dcopy(3,&coords[3*3],1,&coords_edges[2* 3*2+3],1);
  cblas_dcopy(3,&coords[1*3],1,&coords_edges[3* 3*2+0],1); 
  cblas_dcopy(3,&coords[2*3],1,&coords_edges[3* 3*2+3],1); 
  cblas_dcopy(3,&coords[1*3],1,&coords_edges[4* 3*2+0],1);
  cblas_dcopy(3,&coords[3*3],1,&coords_edges[4* 3*2+3],1);
  cblas_dcopy(3,&coords[2*3],1,&coords_edges[5* 3*2+0],1);
  cblas_dcopy(3,&coords[3*3],1,&coords_edges[5* 3*2+3],1);
  PetscFunctionReturn(0);
}
PetscErrorCode FEMethod_ComplexForLazy::GetFint() {
  PetscFunctionBegin;
  try {
  switch (fe_ent_ptr->get_ent_type()) {
    case MBTET: {
      int ee = 0;
      for(;ee<6;ee++) {
  	diff_edgeNinvJac[ee] = &(diffH1edgeNinvJac[ee])[0]; 
      }
      int ff = 0;
      for(;ff<4;ff++) {
        diff_faceNinvJac[ff] = &(diffH1faceNinvJac[ff])[0];
      }
      diff_volumeNinvJac = &diffH1elemNinvJac[0];
      int g_dim = get_dim_gNTET();
      if(!dofs_x_edge_data.empty()) {
	if(dofs_x_edge_data.size()!=6) SETERRQ(PETSC_COMM_SELF,1,"data vectors are not set");
	ee = 0;
	for(;ee<6;ee++) {
	  if(dofs_x_edge_data[ee].size()!=0) {
	    Fint_h_edge_data[ee].resize(dofs_x_edge_data[ee].size());
	    Fint_h_edge[ee] = &Fint_h_edge_data[ee].data()[0];
	  } else {
	    Fint_h_edge[ee] = NULL;
	  }
	}
      }
      if(!dofs_x_face_data.empty()) {
	if(dofs_x_face_data.size()!=4) SETERRQ(PETSC_COMM_SELF,1,"data vectors are not set");
	ff = 0;
	for(;ff<4;ff++) {
	  if(dofs_x_face_data[ff].size()!=0) {
	    Fint_h_face_data[ff].resize(dofs_x_face_data[ff].size());
	    Fint_h_face[ff] = &Fint_h_face_data[ff].data()[0];
	  } else {
	    Fint_h_face[ff] = NULL;
	  }
	}
      }
      if(dofs_x_volume.size()!=0) {
	  assert(RowGlobSpatial[i_volume].size() == (unsigned int)3*NBVOLUME_H1(order_volume));
	  Fint_h_volume.resize(dofs_x_volume.size());
      }
      ierr = GetDofs_X_FromElementData(); CHKERRQ(ierr);
      unsigned int sub_analysis_type[3] = { spatail_analysis, material_analysis, mesh_quality_analysis };
      double _lambda,_mu;
      if( (spatail_analysis|material_analysis)&type_of_analysis ) {
	ierr = GetMatParameters(&_lambda,&_mu,ptr_matctx); CHKERRQ(ierr);
      }
      for(int ss = 0;ss<3;ss++) {
	switch(sub_analysis_type[ss]&type_of_analysis) {
	  case spatail_analysis: {
	    ierr = Fint_Hh_hierarchical(&order_edges[0],&order_faces[0],order_volume,V,_lambda,_mu,ptr_matctx, 
	      &diffNTETinvJac[0],&diff_edgeNinvJac[0],&diff_faceNinvJac[0],&diff_volumeNinvJac[0], 
	      &dofs_X.data()[0],&*dofs_x.data().begin(),NULL,NULL,
	      &dofs_x_edge[0],&dofs_x_face[0],&*dofs_x_volume.data().begin(), 
	      NULL,&*Fint_h.data().begin(),Fint_h_edge,Fint_h_face,&*Fint_h_volume.data().begin(),
	      NULL,NULL,NULL,NULL,NULL,
	      g_dim,g_TET_W); CHKERRQ(ierr);
	  }
	  break;
	  case material_analysis: {
	    ierr = Fint_Hh_hierarchical(&order_edges[0],&order_faces[0],order_volume,V,_lambda,_mu,ptr_matctx, 
	      &diffNTETinvJac[0],&diff_edgeNinvJac[0],&diff_faceNinvJac[0],&diff_volumeNinvJac[0], 
	      &dofs_X.data()[0],&*dofs_x.data().begin(),NULL,NULL,
	      &dofs_x_edge[0],&dofs_x_face[0],&*dofs_x_volume.data().begin(), 
	      &*Fint_H.data().begin(),NULL,NULL,NULL,NULL,
	      NULL,NULL,NULL,NULL,NULL,
	      g_dim,g_TET_W); CHKERRQ(ierr);
	  }
	  break;
	  case mesh_quality_analysis: {
	    if(!flg_alpha2) SETERRQ(PETSC_COMM_SELF,1,"-my_alpha2 is not set");
	    if(!flg_gamma) SETERRQ(PETSC_COMM_SELF,1,"-my_gamma is not set");
	    double coords_edges[2*3*6]; 
	    double alpha2_array[4] = { alpha2,alpha2,alpha2,alpha2 }; 
	    ierr = get_edges_from_elem_coords(&coords[0],coords_edges); CHKERRQ(ierr);
	    ierr = quality_volume_length_F(V,alpha2_array,gamma,&diffNTETinvJac[0],coords_edges,
	      &dofs_X.data()[0],NULL,NULL,NULL,quality0,quality,b,
	      &*Fint_H.data().begin(),NULL); CHKERRQ(ierr);
	  }
	  break;
	  default: 
	    continue;
      }
    }
    } break;
    default:
      SETERRQ(PETSC_COMM_SELF,1,"not implemented");
  }
  } catch (const std::exception& ex) {
    ostringstream ss;
    ss << "throw in method: " << ex.what() << endl;
    SETERRQ(PETSC_COMM_SELF,1,ss.str().c_str());
  } 
  PetscFunctionReturn(0);
}
PetscErrorCode FEMethod_ComplexForLazy::GetFaceIndicesAndData(EntityHandle face) {
  PetscFunctionBegin;
  typedef FENumeredDofMoFEMEntity_multiIndex::index<Composite_Name_And_Ent_mi_tag>::type::iterator dofs_iterator;
  try {
  FaceNodeIndices.resize(9);
  fill(FaceNodeIndices.begin(),FaceNodeIndices.end(),-1);
  FaceNodeData.resize(9);
  const EntityHandle* conn_face; 
  int num_nodes; 
  rval = moab.get_connectivity(face,conn_face,num_nodes,true); CHKERR_PETSC(rval);
  int nn = 0,dd = 0;
  for(;nn<3;nn++) {
    dofs_iterator niit,hi_niit;
    dofs_iterator col_niit,hi_col_niit;
    niit = row_multiIndex->get<Composite_Name_And_Ent_mi_tag>().lower_bound(boost::make_tuple(spatial_field_name,conn_face[nn]));
    hi_niit = row_multiIndex->get<Composite_Name_And_Ent_mi_tag>().upper_bound(boost::make_tuple(spatial_field_name,conn_face[nn]));
    col_niit = col_multiIndex->get<Composite_Name_And_Ent_mi_tag>().lower_bound(boost::make_tuple(spatial_field_name,conn_face[nn]));
    hi_col_niit = col_multiIndex->get<Composite_Name_And_Ent_mi_tag>().upper_bound(boost::make_tuple(spatial_field_name,conn_face[nn]));
    for(;niit!=hi_niit;niit++,col_niit++,dd++) {
      if(col_niit->get_petsc_gloabl_dof_idx() != niit->get_petsc_gloabl_dof_idx()) {
	SETERRQ(PETSC_COMM_SELF,1,"this implementation is not working for diffrent ordering of columns and rows"); 
      }
      FaceNodeIndices[nn*niit->get_max_rank()+niit->get_dof_rank()] = niit->get_petsc_gloabl_dof_idx();
      FaceNodeData[nn*niit->get_max_rank()+niit->get_dof_rank()] = niit->get_FieldData();
    }
  }
  if(dd != 9) {
    SETERRQ(PETSC_COMM_SELF,1,"face is not adjacent to this TET"); 
  }
  dofs_iterator fiit,hi_fiit;
  fiit = row_multiIndex->get<Composite_Name_And_Ent_mi_tag>().lower_bound(boost::make_tuple(spatial_field_name,face));
  hi_fiit = row_multiIndex->get<Composite_Name_And_Ent_mi_tag>().upper_bound(boost::make_tuple(spatial_field_name,face));
  if(fiit!=hi_fiit) {
    FaceIndices.resize(distance(fiit,hi_fiit));
    FaceData.resize(distance(fiit,hi_fiit));
    face_order = fiit->get_max_order();
    assert((unsigned int)3*NBFACE_H1(face_order)==distance(fiit,hi_fiit));
    if(NBFACE_H1(face_order)>0) {
      dd = 0;
      for(dofs_iterator fiiit = fiit;fiiit!=hi_fiit;fiiit++,dd++) {
	FaceIndices[fiiit->get_EntDofIdx()] = fiiit->get_petsc_gloabl_dof_idx();
	FaceData[fiiit->get_EntDofIdx()] = fiiit->get_FieldData();
      }
    }
    N_face.resize(g_TRI_dim*NBFACE_H1(face_order));
    diffN_face.resize(2*g_TRI_dim*NBFACE_H1(face_order));
    int face_nodes[] = { 0,1,2 };
    ierr = H1_FaceShapeFunctions_MBTRI(face_nodes,face_order,&g_NTRI[0],&diffNTRI[0],&N_face[0],&diffN_face[0],g_TRI_dim); CHKERRQ(ierr);
  } else {
    face_order = 0;
  }
  FaceEdgeIndices_data.resize(3);
  FaceEdgeData_data.resize(3);
  FaceEdgeSense.resize(3);
  FaceEdgeOrder.resize(3);
  N_edge_data.resize(3);
  diffN_edge_data.resize(3);
  int ee = 0;
  for(;ee<3;ee++) {
    EntityHandle edge;
    rval = moab.side_element(face,1,ee,edge); CHKERR_PETSC(rval);
    int side_number,offset;
    rval = moab.side_number(face,edge,side_number,FaceEdgeSense[ee],offset); CHKERR_PETSC(rval);
    dofs_iterator eiit,hi_eiit;
    eiit = row_multiIndex->get<Composite_Name_And_Ent_mi_tag>().lower_bound(boost::make_tuple(spatial_field_name,edge));
    hi_eiit = row_multiIndex->get<Composite_Name_And_Ent_mi_tag>().upper_bound(boost::make_tuple(spatial_field_name,edge));
    if(eiit!=hi_eiit) {
      FaceEdgeOrder[ee] = eiit->get_max_order();
      if(NBEDGE_H1(FaceEdgeOrder[ee])>0) {
	assert(3*NBEDGE_H1(FaceEdgeOrder[ee]) == distance(eiit,hi_eiit));
	FaceEdgeIndices_data[ee].resize(distance(eiit,hi_eiit));
	FaceEdgeData_data[ee].resize(distance(eiit,hi_eiit));
	for(;eiit!=hi_eiit;eiit++) {
	  FaceEdgeIndices_data[ee][eiit->get_EntDofIdx()] = eiit->get_petsc_gloabl_dof_idx();
	  FaceEdgeData_data[ee][eiit->get_EntDofIdx()] = eiit->get_FieldData();
	}
	EdgeData[ee] = &*(FaceEdgeData_data[ee].data().begin());
	N_edge_data[ee].resize(g_TRI_dim*NBEDGE_H1(FaceEdgeOrder[ee]));
	diffN_edge_data[ee].resize(2*g_TRI_dim*NBEDGE_H1(FaceEdgeOrder[ee]));
	N_edge[ee] = &(N_edge_data[ee][0]);
	diffN_edge[ee] = &(diffN_edge_data[ee][0]);
      }
    } else {
      FaceEdgeOrder[ee] = 0;
      EdgeData[ee] = NULL;
      N_edge[ee] = NULL;
      diffN_edge[ee] = NULL;
    }
  }
  ierr = H1_EdgeShapeFunctions_MBTRI(&FaceEdgeSense[0],&FaceEdgeOrder[0],&g_NTRI[0],&diffNTRI[0],N_edge,diffN_edge,g_TRI_dim); CHKERRQ(ierr);
  GetFaceIndicesAndData_face = face;
  } catch (const std::exception& ex) {
    ostringstream ss;
    ss << "throw in method: " << ex.what() << endl;
    SETERRQ(PETSC_COMM_SELF,1,ss.str().c_str());
  } 
  PetscFunctionReturn(0);
}
PetscErrorCode FEMethod_ComplexForLazy::GetFExt(EntityHandle face,double *t,double *t_edge[],double *t_face) {
  PetscFunctionBegin;
  try {
  if(GetFaceIndicesAndData_face!=face) SETERRQ(PETSC_COMM_SELF,1,"run GetFaceIndicesAndData(face) before call of this function");
  FExt.resize(9);
  FExt_edge_data.resize(3);
  int ee = 0;
  for(;ee<3;ee++) {
    if(FaceEdgeIndices_data[ee].size()>0) {
      if((unsigned int)3*NBEDGE_H1(FaceEdgeOrder[ee])!=FaceEdgeIndices_data[ee].size()) SETERRQ(PETSC_COMM_SELF,1,"data inconsistency"); 
      FExt_edge_data[ee].resize(FaceEdgeIndices_data[ee].size());
      FExt_edge[ee] = &*(FExt_edge_data[ee].data().begin());
    } else {
      if(NBEDGE_H1(FaceEdgeOrder[ee])!=0) SETERRQ(PETSC_COMM_SELF,1,"data inconsistency"); 
      FExt_edge[ee] = NULL;
    }
  }
  FExt_face.resize(FaceIndices.size());
  /*const EntityHandle* conn_face;
  ublas::vector<double> coords_face;
  coords_face.resize(9);
  int num_nodes;
  rval = moab.get_connectivity(face,conn_face,num_nodes,true); CHKERR_PETSC(rval);
  rval = moab.get_coords(conn_face,num_nodes,&*coords_face.begin()); CHKERR_PETSC(rval);*/
  switch(type_of_forces) {
    case conservative:
      ierr = Fext_h_hierarchical(
	face_order,&FaceEdgeOrder[0],//2
	&g_NTRI[0],&N_face[0],N_edge,&diffNTRI[0],&diffN_face[0],diffN_edge,//8
	t,t_edge,t_face,//11
	&*FaceNodeData.data().begin(),EdgeData,&*FaceData.data().begin(),//14
	//&*coords_face.data().begin(),NULL,NULL,//14
	NULL,NULL,NULL,//17
	&*FExt.data().begin(),FExt_edge,&*FExt_face.data().begin(),//20
	NULL,NULL,NULL,//23
	g_TRI_dim,g_TRI_W); CHKERRQ(ierr);
      break;
    case nonconservative:
      SETERRQ(PETSC_COMM_SELF,1,"not implemented"); 
      break;
  }
  } catch (const std::exception& ex) {
      ostringstream ss;
      ss << "throw in method: " << ex.what() << endl;
      SETERRQ(PETSC_COMM_SELF,1,ss.str().c_str());
  } 
  PetscFunctionReturn(0);
}
PetscErrorCode FEMethod_ComplexForLazy::GetTangentExt(EntityHandle face,double *t,double *t_edge[],double *t_face) {
  PetscFunctionBegin;
  try {
  if(GetFaceIndicesAndData_face!=face) SETERRQ(PETSC_COMM_SELF,1,"run GetFaceIndicesAndData(face) before call of this function");
  switch(type_of_forces) {
    case conservative: {
      const EntityHandle* conn_face;
      ublas::vector<double> coords_face;
      coords_face.resize(9);
      int num_nodes;
      rval = moab.get_connectivity(face,conn_face,num_nodes,true); CHKERR_PETSC(rval);
      rval = moab.get_coords(conn_face,num_nodes,&*coords_face.begin()); CHKERR_PETSC(rval);
      double center[3]; 
      tricircumcenter3d_tp(&coords_face.data()[0],&coords_face.data()[3],&coords_face.data()[6],center,NULL,NULL);
      cblas_daxpy(3,-1,&coords_face.data()[0],1,center,1);
      double r = cblas_dnrm2(3,center,1);
      KExt_hh.resize(9,9);
      KExt_edgeh_data.resize(3);
      int ee = 0;
      for(;ee<3;ee++) {
        assert((unsigned int)3*NBEDGE_H1(FaceEdgeOrder[ee]) == FaceEdgeIndices_data[ee].size());
        KExt_edgeh_data[ee].resize(FaceEdgeIndices_data[ee].size(),9);
        KExt_edgeh[ee] = &*KExt_edgeh_data[ee].data().begin();
      }
      KExt_faceh.resize(FaceIndices.size(),9);
      ierr = KExt_hh_hierarchical(r*eps,face_order,&FaceEdgeOrder[0],
          &g_NTRI[0],&N_face[0],N_edge,&diffNTRI[0],&diffN_face[0],diffN_edge,
          t,t_edge,t_face,&*FaceNodeData.data().begin(),EdgeData,&*FaceData.data().begin(),
          &*KExt_hh.data().begin(),KExt_edgeh,&*KExt_faceh.data().begin(),g_TRI_dim,g_TRI_W); CHKERRQ(ierr);
      //
      KExt_hedge_data.resize(3);
      KExt_edgeedge_data.resize(3,3);
      KExt_faceedge_data.resize(3);
      for(ee = 0;ee<3;ee++) {
        KExt_hedge_data[ee].resize(9,FaceEdgeIndices_data[ee].size());
        KExt_hedge[ee] = &*KExt_hedge_data[ee].data().begin();
        for(int eee = 0;eee<3;eee++) {
          KExt_edgeedge_data(ee,eee).resize(FaceEdgeIndices_data[ee].size(),FaceEdgeIndices_data[eee].size());
          KExt_edgeedge[ee][eee] = &*KExt_edgeedge_data(ee,eee).data().begin();
        }
        KExt_faceedge_data[ee].resize(FaceIndices.size(),FaceEdgeIndices_data[ee].size());
        KExt_faceedge[ee] = &*KExt_faceedge_data[ee].data().begin();
      }
      ierr = KExt_hh_hierarchical_edge(r*eps,face_order,&FaceEdgeOrder[0],
          &g_NTRI[0],&N_face[0],N_edge,&diffNTRI[0],&diffN_face[0],diffN_edge,
          t,t_edge,t_face,&*FaceNodeData.data().begin(),EdgeData,&*FaceData.data().begin(),
          KExt_hedge,KExt_edgeedge,KExt_faceedge,
          g_TRI_dim,g_TRI_W); CHKERRQ(ierr);
      //
      KExt_hface.resize(9,FaceIndices.size());
      KExt_edgeface_data.resize(3);
      for(ee = 0;ee<3;ee++) {
        KExt_edgeface_data[ee].resize(FaceEdgeIndices_data[ee].size(),FaceIndices.size());
        KExt_edgeface[ee] = &*KExt_edgeface_data[ee].data().begin();
      }
      assert(FaceIndices.size() == (unsigned int)3*NBFACE_H1(face_order));
      KExt_faceface.resize(FaceIndices.size(),FaceIndices.size());
      ierr = KExt_hh_hierarchical_face(r*eps,face_order,&FaceEdgeOrder[0],
          &g_NTRI[0],&N_face[0],N_edge,&diffNTRI[0],&diffN_face[0],diffN_edge,
          t,t_edge,t_face,&*FaceNodeData.data().begin(),EdgeData,&*FaceData.data().begin(),
          &*KExt_hface.data().begin(),KExt_edgeface,&*KExt_faceface.data().begin(),
          g_TRI_dim,g_TRI_W); CHKERRQ(ierr); 
      }
      break;
    case nonconservative:
      SETERRQ(PETSC_COMM_SELF,1,"not implemented"); 
      break;
    }
  } catch (const std::exception& ex) {
      ostringstream ss;
      ss << "throw in method: " << ex.what() << endl;
      SETERRQ(PETSC_COMM_SELF,1,ss.str().c_str());
  } 
  PetscFunctionReturn(0);
}
PetscErrorCode FEMethod_ComplexForLazy::GetFaceIndicesAndData_Material(EntityHandle face) {
  PetscFunctionBegin;
  typedef FENumeredDofMoFEMEntity_multiIndex::index<Composite_Name_And_Ent_mi_tag>::type::iterator dofs_iterator;
  typedef FEDofMoFEMEntity_multiIndex::index<Composite_Name_And_Ent_mi_tag>::type::iterator data_dofs_iterator;
  try {
  FaceNodeIndices_Material.resize(9);
  fill(FaceNodeIndices_Material.begin(),FaceNodeIndices_Material.end(),-1);
  FaceNodeData_Material.resize(9);
  const EntityHandle* conn_face; 
  int num_nodes; 
  rval = moab.get_connectivity(face,conn_face,num_nodes,true); CHKERR_PETSC(rval);
  int nn = 0,dd = 0;
  for(;nn<3;nn++) {
    dofs_iterator niit,hi_niit;
    dofs_iterator col_niit,hi_col_niit;
    niit = row_multiIndex->get<Composite_Name_And_Ent_mi_tag>().lower_bound(boost::make_tuple(material_field_name,conn_face[nn]));
    hi_niit = row_multiIndex->get<Composite_Name_And_Ent_mi_tag>().upper_bound(boost::make_tuple(material_field_name,conn_face[nn]));
    col_niit = col_multiIndex->get<Composite_Name_And_Ent_mi_tag>().lower_bound(boost::make_tuple(material_field_name,conn_face[nn]));
    hi_col_niit = col_multiIndex->get<Composite_Name_And_Ent_mi_tag>().upper_bound(boost::make_tuple(material_field_name,conn_face[nn]));
    for(;niit!=hi_niit;niit++,col_niit++,dd++) {
      assert(col_niit->get_petsc_gloabl_dof_idx() == niit->get_petsc_gloabl_dof_idx());
      FaceNodeIndices_Material[nn*niit->get_max_rank()+niit->get_EntDofIdx()] = niit->get_petsc_gloabl_dof_idx();
      FaceNodeData_Material[nn*niit->get_max_rank()+niit->get_EntDofIdx()] = niit->get_FieldData();
    }
  }
  if(dd != 9) {
    SETERRQ(PETSC_COMM_SELF,1,"face is not adjacent to this TET"); 
  }
  FaceNodeData.resize(9);
  nn = dd = 0;
  for(;nn<3;nn++) {
    data_dofs_iterator niit,hi_niit;
    niit = data_multiIndex->get<Composite_Name_And_Ent_mi_tag>().lower_bound(boost::make_tuple(spatial_field_name,conn_face[nn]));
    hi_niit = data_multiIndex->get<Composite_Name_And_Ent_mi_tag>().upper_bound(boost::make_tuple(spatial_field_name,conn_face[nn]));
    for(;niit!=hi_niit;niit++,dd++) {
      FaceNodeData[nn*niit->get_max_rank()+niit->get_EntDofIdx()] = niit->get_FieldData();
    }
  }
  if(dd != 9) {
    SETERRQ(PETSC_COMM_SELF,1,"face is not adjacent to this TET"); 
  }
  data_dofs_iterator fiit,hi_fiit;
  fiit = data_multiIndex->get<Composite_Name_And_Ent_mi_tag>().lower_bound(boost::make_tuple(spatial_field_name,face));
  hi_fiit = data_multiIndex->get<Composite_Name_And_Ent_mi_tag>().upper_bound(boost::make_tuple(spatial_field_name,face));
  if(fiit!=hi_fiit) {
    FaceData.resize(distance(fiit,hi_fiit));
    face_order = fiit->get_max_order();
    assert((unsigned int)3*NBFACE_H1(face_order)==distance(fiit,hi_fiit));
    if(NBFACE_H1(face_order)>0) {
      for(data_dofs_iterator fiiit = fiit;fiiit!=hi_fiit;fiiit++) {
	FaceData[fiiit->get_EntDofIdx()] = fiiit->get_FieldData();
      }
    }
    N_face.resize(g_TRI_dim*NBFACE_H1(face_order));
    diffN_face.resize(2*g_TRI_dim*NBFACE_H1(face_order));
    int face_nodes[] = { 0,1,2 };
    ierr = H1_FaceShapeFunctions_MBTRI(face_nodes,face_order,&g_NTRI[0],&diffNTRI[0],&N_face[0],&diffN_face[0],g_TRI_dim); CHKERRQ(ierr);
  } else {
    face_order = 0;
  }
  FaceEdgeData_data.resize(3);
  FaceEdgeSense.resize(3);
  FaceEdgeOrder.resize(3);
  N_edge_data.resize(3);
  diffN_edge_data.resize(3);
  int ee = 0;
  for(;ee<3;ee++) {
    EntityHandle edge;
    rval = moab.side_element(face,1,ee,edge); CHKERR_PETSC(rval);
    int side_number,offset;
    rval = moab.side_number(face,edge,side_number,FaceEdgeSense[ee],offset); CHKERR_PETSC(rval);
    data_dofs_iterator eiit,hi_eiit;
    eiit = data_multiIndex->get<Composite_Name_And_Ent_mi_tag>().lower_bound(boost::make_tuple(spatial_field_name,edge));
    hi_eiit = data_multiIndex->get<Composite_Name_And_Ent_mi_tag>().upper_bound(boost::make_tuple(spatial_field_name,edge));
    if(eiit!=hi_eiit) {
      FaceEdgeOrder[ee] = eiit->get_max_order();
      if(NBEDGE_H1(FaceEdgeOrder[ee])>0) {
	assert(3*NBEDGE_H1(FaceEdgeOrder[ee]) == distance(eiit,hi_eiit));
	FaceEdgeData_data[ee].resize(distance(eiit,hi_eiit));
	for(;eiit!=hi_eiit;eiit++) {
	  FaceEdgeData_data[ee][eiit->get_EntDofIdx()] = eiit->get_FieldData();
	}
	EdgeData[ee] = &*(FaceEdgeData_data[ee].data().begin());
	N_edge_data[ee].resize(g_TRI_dim*NBEDGE_H1(FaceEdgeOrder[ee]));
	diffN_edge_data[ee].resize(2*g_TRI_dim*NBEDGE_H1(FaceEdgeOrder[ee]));
	N_edge[ee] = &(N_edge_data[ee][0]);
	diffN_edge[ee] = &(diffN_edge_data[ee][0]);
      }
    } else {
      FaceEdgeOrder[ee] = 0;
      EdgeData[ee] = NULL;
    }
  }
  ierr = H1_EdgeShapeFunctions_MBTRI(&FaceEdgeSense[0],&FaceEdgeOrder[0],&g_NTRI[0],&diffNTRI[0],N_edge,diffN_edge,g_TRI_dim); CHKERRQ(ierr);
  GetFaceIndicesAndData_face = face;
  } catch (const std::exception& ex) {
    ostringstream ss;
    ss << "throw in method: " << ex.what() << endl;
    SETERRQ(PETSC_COMM_SELF,1,ss.str().c_str());
  } 
  PetscFunctionReturn(0);
}
PetscErrorCode FEMethod_ComplexForLazy::GetFExt_Material(EntityHandle face,double *t,double *t_edge[],double *t_face) {
  PetscFunctionBegin;
  try {
  if(GetFaceIndicesAndData_face!=face) SETERRQ(PETSC_COMM_SELF,1,"run GetFaceIndicesAndData_Material(face) before call of this function");
  FExt_Material.resize(9);
  switch(type_of_forces) {
    case conservative:
      ierr = Fext_H(
	face_order,&FaceEdgeOrder[0],//2
	&g_NTRI[0],&N_face[0],N_edge,&diffNTRI[0],&diffN_face[0],diffN_edge,//8
	t,t_edge,t_face,//11
	&*FaceNodeData_Material.data().begin(),NULL,
	&*FaceNodeData.data().begin(),EdgeData,&*FaceData.data().begin(),
	NULL,NULL,NULL,
	&*FExt_Material.begin(),NULL,g_TRI_dim,g_TRI_W); CHKERRQ(ierr);
      break;
    case nonconservative:
      SETERRQ(PETSC_COMM_SELF,1,"not implemented"); 
      break;
  }
  FExt_Material *= -1;
  } catch (const std::exception& ex) {
      ostringstream ss;
      ss << "throw in method: " << ex.what() << endl;
      SETERRQ(PETSC_COMM_SELF,1,ss.str().c_str());
  } 
  PetscFunctionReturn(0);
}
PetscErrorCode FEMethod_ComplexForLazy::GetTangentExt_Material(EntityHandle face,double *t,double *t_edge[],double *t_face) {
  PetscFunctionBegin;
  try {
    switch(type_of_forces) {
      case conservative: {
        if(GetFaceIndicesAndData_face!=face) SETERRQ(PETSC_COMM_SELF,1,"run GetFaceIndicesAndData_Material(face) before call of this function");
        const EntityHandle* conn_face;
        ublas::vector<double> coords_face;
        coords_face.resize(9);
        int num_nodes;
        rval = moab.get_connectivity(face,conn_face,num_nodes,true); CHKERR_PETSC(rval);
        rval = moab.get_coords(conn_face,num_nodes,&*coords_face.begin()); CHKERR_PETSC(rval);
	double center[3]; 
	tricircumcenter3d_tp(&coords_face.data()[0],&coords_face.data()[3],&coords_face.data()[6],center,NULL,NULL);
	cblas_daxpy(3,-1,&coords_face.data()[0],1,center,1);
	double r = cblas_dnrm2(3,center,1);
        KExt_HH_Material.resize(9,9);
        ierr = KExt_HH(r*eps,face_order,&FaceEdgeOrder[0],
          &g_NTRI[0],&N_face[0],N_edge,&diffNTRI[0],&diffN_face[0],diffN_edge,
          t,t_edge,t_face,
          &*FaceNodeData_Material.data().begin(),
          &*FaceNodeData.data().begin(),EdgeData,&*FaceData.data().begin(),
          &*KExt_HH_Material.data().begin(),g_TRI_dim,g_TRI_W); CHKERRQ(ierr);
      }
      break;
    case nonconservative:
      SETERRQ(PETSC_COMM_SELF,1,"not implemented"); 
      break;
    }
  } catch (const std::exception& ex) {
      ostringstream ss;
      ss << "throw in method: " << ex.what() << endl;
      SETERRQ(PETSC_COMM_SELF,1,ss.str().c_str());
  } 
  PetscFunctionReturn(0);
}


}

