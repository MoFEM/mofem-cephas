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

#include "moabFEMethod_ComplexForLazy.hpp"
#include "FEM.h"
#include "complex_for_lazy.h"

extern "C" {
void tetcircumcenter_tp(double a[3],double b[3],double c[3], double d[3],
  double circumcenter[3],double *xi,double *eta,double *zeta);
void tricircumcenter3d_tp(double a[3],double b[3],double c[3],
  double circumcenter[3],double *xi,double *eta);
}

namespace MoFEM {

FEMethod_ComplexForLazy::FEMethod_ComplexForLazy(Interface& _moab,analysis _type,
    double _lambda,double _mu, int _verbose): 
    FEMethod_UpLevelStudent(_moab,_verbose), type_of_analysis(_type), lambda(_lambda),mu(_mu), eps(1e-12) {
  pcomm = ParallelComm::get_pcomm(&moab,MYPCOMM_INDEX);
  order_edges.resize(6);
  order_faces.resize(4);
  edgeNinvJac.resize(6);
  faceNinvJac.resize(4);
  diff_edgeNinvJac.resize(6);
  diff_faceNinvJac.resize(4);
  //Tangent_hh_hierachical
  Kedgeh_data.resize(6);
  Kfaceh_data.resize(4);
  //Tangent_hh_hierachical_edge
  Khedge_data.resize(6);
  Khh_volumeedge_data.resize(6);
  Khh_edgeedge_data.resize(6,6);
  Khh_faceedge_data.resize(4,6);
  //Tangent_hh_hierachical_face
  Khface_data.resize(6);
  Khh_volumeface_data.resize(6);
  Khh_faceface_data.resize(4,4);
  Khh_edgeface_data.resize(6,4);
  //
  dofs_x.resize(12);
  dofs_x_edge_data.resize(6);
  dofs_x_face_data.resize(4);
  dofs_x_edge.resize(6);
  dofs_x_face.resize(4);
  //
  Khh_edgevolume_data.resize(6);
  Khh_facevolume_data.resize(4);
  //
  Fblock_x.resize(12);
  Fint_h_edge_data.resize(6);
  Fint_h_face_data.resize(4);
  //
  g_NTET.resize(4*45);
  ShapeMBTET(&g_NTET[0],G_TET_X45,G_TET_Y45,G_TET_Z45,45);
  g_TET_W = G_TET_W45;
  //
  g_NTRI.resize(3*7);
  ShapeMBTRI(&g_NTRI[0],G_TRI_X7,G_TRI_Y7,7); 
  g_TRI_dim = 7;
  g_TRI_W = G_TRI_W7;
}
PetscErrorCode FEMethod_ComplexForLazy::OpComplexForLazyStart() {
  PetscFunctionBegin;
  ierr = OpStudentStart_TET(g_NTET); CHKERRQ(ierr);
  PetscFunctionReturn(0);
}
PetscErrorCode FEMethod_ComplexForLazy::GetIndices() {
  PetscFunctionBegin;
  switch (fe_ent_ptr->get_ent_type()) {
    case MBTET: {
      RowGlob.resize(1+6+4+1);
      ColGlob.resize(1+6+4+1);
      rowNMatrices.resize(1+6+4+1);
      rowDiffNMatrices.resize(1+6+4+1);
      rowBMatrices.resize(1+6+4+1);
      colNMatrices.resize(1+6+4+1);
      colDiffNMatrices.resize(1+6+4+1);
      colBMatrices.resize(1+6+4+1);
      try {
      //nodes
      ierr = GetRowIndices("SPATIAL_POSITION",RowGlob[0]); CHKERRQ(ierr);
      ierr = GetColIndices("SPATIAL_POSITION",ColGlob[0]); CHKERRQ(ierr);
      //row
      ierr = GetGaussRowNMatrix("SPATIAL_POSITION",rowNMatrices[0]); CHKERRQ(ierr);
      ierr = GetGaussRowDiffNMatrix("SPATIAL_POSITION",rowDiffNMatrices[0]); CHKERRQ(ierr);
      ierr = MakeBMatrix3D("SPATIAL_POSITION",rowDiffNMatrices[0],rowBMatrices[0]);  CHKERRQ(ierr);
      //col
      ierr = GetGaussColNMatrix("SPATIAL_POSITION",colNMatrices[0]); CHKERRQ(ierr);
      ierr = GetGaussColDiffNMatrix("SPATIAL_POSITION",colDiffNMatrices[0]); CHKERRQ(ierr);
      ierr = MakeBMatrix3D("SPATIAL_POSITION",colDiffNMatrices[0],colBMatrices[0]);  CHKERRQ(ierr);
      //edges
      int ee = 0;
      for(;ee<6;ee++) { //edges matrices
	ierr = GetRowIndices("SPATIAL_POSITION",MBEDGE,RowGlob[1+ee],ee); CHKERRQ(ierr);
	if(RowGlob[1+ee].size()!=0) {
	  ierr = GetGaussRowNMatrix("SPATIAL_POSITION",MBEDGE,rowNMatrices[1+ee],ee); CHKERRQ(ierr);
	  ierr = GetGaussRowDiffNMatrix("SPATIAL_POSITION",MBEDGE,rowDiffNMatrices[1+ee],ee); CHKERRQ(ierr);
	  ierr = MakeBMatrix3D("SPATIAL_POSITION",rowDiffNMatrices[1+ee],rowBMatrices[1+ee]);  CHKERRQ(ierr);
	}
	ierr = GetColIndices("SPATIAL_POSITION",MBEDGE,ColGlob[1+ee],ee); CHKERRQ(ierr);
	if(ColGlob[1+ee].size()>0) {
	  ierr = GetGaussColNMatrix("SPATIAL_POSITION",MBEDGE,colNMatrices[1+ee],ee); CHKERRQ(ierr);
	  ierr = GetGaussColDiffNMatrix("SPATIAL_POSITION",MBEDGE,colDiffNMatrices[1+ee],ee); CHKERRQ(ierr);
	  ierr = MakeBMatrix3D("SPATIAL_POSITION",colDiffNMatrices[1+ee],colBMatrices[1+ee]);  CHKERRQ(ierr);
	}
      }
      assert(ee == 6);
      //faces
      int ff = 0;
      for(;ff<4;ff++) { //faces matrices
	ierr = GetRowIndices("SPATIAL_POSITION",MBTRI,RowGlob[1+ee+ff],ff); CHKERRQ(ierr);
	if(RowGlob[1+ee+ff].size()!=0) {
	  ierr = GetGaussRowNMatrix("SPATIAL_POSITION",MBTRI,rowNMatrices[1+ee+ff],ff); CHKERRQ(ierr);
	  ierr = GetGaussRowDiffNMatrix("SPATIAL_POSITION",MBTRI,rowDiffNMatrices[1+ee+ff],ff); CHKERRQ(ierr);
	  ierr = MakeBMatrix3D("SPATIAL_POSITION",rowDiffNMatrices[1+ee+ff],rowBMatrices[1+ee+ff]);  CHKERRQ(ierr);
	}
	ierr = GetColIndices("SPATIAL_POSITION",MBTRI,ColGlob[1+ee+ff],ff); CHKERRQ(ierr);
	if(ColGlob[1+ee+ff].size()!=0) {
	  ierr = GetGaussColNMatrix("SPATIAL_POSITION",MBTRI,colNMatrices[1+ee+ff],ff); CHKERRQ(ierr);
	  ierr = GetGaussColDiffNMatrix("SPATIAL_POSITION",MBTRI,colDiffNMatrices[1+ee+ff],ff); CHKERRQ(ierr);
	  ierr = MakeBMatrix3D("SPATIAL_POSITION",colDiffNMatrices[1+ee+ff],colBMatrices[1+ee+ff]);  CHKERRQ(ierr);
	}
      }
      assert(ff == 4);
      //volumes
      ierr = GetRowIndices("SPATIAL_POSITION",MBTET,RowGlob[1+ee+ff]); CHKERRQ(ierr);
      if(RowGlob[1+ee+ff].size()!=0) {
	ierr = GetGaussRowNMatrix("SPATIAL_POSITION",MBTET,rowNMatrices[1+ee+ff]); CHKERRQ(ierr);
	ierr = GetGaussRowDiffNMatrix("SPATIAL_POSITION",MBTET,rowDiffNMatrices[1+ee+ff]); CHKERRQ(ierr);
	ierr = MakeBMatrix3D("SPATIAL_POSITION",rowDiffNMatrices[1+ee+ff],rowBMatrices[1+ee+ff]);  CHKERRQ(ierr);
      }
      ierr = GetColIndices("SPATIAL_POSITION",MBTET,ColGlob[1+ee+ff]); CHKERRQ(ierr);
      if(ColGlob[1+ee+ff].size()!=0) {
	ierr = GetGaussColNMatrix("SPATIAL_POSITION",MBTET,colNMatrices[1+ee+ff]); CHKERRQ(ierr);
	ierr = GetGaussColDiffNMatrix("SPATIAL_POSITION",MBTET,colDiffNMatrices[1+ee+ff]); CHKERRQ(ierr);
	ierr = MakeBMatrix3D("SPATIAL_POSITION",colDiffNMatrices[1+ee+ff],colBMatrices[1+ee+ff]);  CHKERRQ(ierr);
      }
      //data edge
      ee = 0;
      for(;ee<6;ee++) {
	FENumeredDofMoFEMEntity_multiIndex::index<Composite_mi_tag>::type::iterator eiit,hi_eiit;
	eiit = row_multiIndex->get<Composite_mi_tag>().lower_bound(boost::make_tuple("SPATIAL_POSITION",MBEDGE,ee));
	if(eiit!=row_multiIndex->get<Composite_mi_tag>().end()) {
	  order_edges[ee] = eiit->get_max_order();
	  hi_eiit = row_multiIndex->get<Composite_mi_tag>().upper_bound(boost::make_tuple("SPATIAL_POSITION",MBEDGE,ee));
	  dofs_x_edge_data[ee].resize(distance(eiit,hi_eiit));
	  dofs_x_edge[ee] = &dofs_x_edge_data[ee].data()[0];
	  assert(dofs_x_edge_data[ee].size() == 3*(unsigned int)NBEDGE_H1(order_edges[ee]));
	  for(int dd = 0;eiit!=hi_eiit;eiit++,dd++) dofs_x_edge_data[ee][dd] = eiit->get_FieldData(); 
	} else {
	  order_edges[ee] = 0;
	}
      }
      //data face
      ff = 0;
      for(;ff<4;ff++) {
	FENumeredDofMoFEMEntity_multiIndex::index<Composite_mi_tag>::type::iterator fiit,hi_fiit;
	fiit = row_multiIndex->get<Composite_mi_tag>().lower_bound(boost::make_tuple("SPATIAL_POSITION",MBTRI,ff));
	if(fiit!=row_multiIndex->get<Composite_mi_tag>().end()) {
	  order_faces[ff] = fiit->get_max_order();
	  hi_fiit = row_multiIndex->get<Composite_mi_tag>().upper_bound(boost::make_tuple("SPATIAL_POSITION",MBTRI,ff));
	  dofs_x_face_data[ff].resize(distance(fiit,hi_fiit));
	  dofs_x_face[ff] = &dofs_x_face_data[ff].data()[0];
	  assert(dofs_x_face_data[ff].size() == 3*(unsigned int)NBFACE_H1(order_faces[ff]));
	  for(int dd = 0;fiit!=hi_fiit;fiit++,dd++) dofs_x_face_data[ff][dd] = fiit->get_FieldData(); 
	} else {
	  order_faces[ee] = 0;
	}
      }
      //data voolume
      FENumeredDofMoFEMEntity_multiIndex::index<Composite_mi_tag2>::type::iterator viit,hi_viit;
      viit = row_multiIndex->get<Composite_mi_tag2>().lower_bound(boost::make_tuple("SPATIAL_POSITION",MBTET));
      order_volume = viit->get_max_order();
      hi_viit = row_multiIndex->get<Composite_mi_tag2>().upper_bound(boost::make_tuple("SPATIAL_POSITION",MBTET));
      dofs_x_volume.resize(distance(viit,hi_viit));
      assert(dofs_x_volume.size() == (unsigned int)NBFACE_H1(order_volume));
      for(int dd = 0;viit!=hi_viit;viit++,dd++) dofs_x_volume[dd] = viit->get_FieldData(); 
      //data nodes
      FENumeredDofMoFEMEntity_multiIndex::index<Composite_mi_tag2>::type::iterator niit,hi_niit;
      niit = row_multiIndex->get<Composite_mi_tag2>().lower_bound(boost::make_tuple("SPATIAL_POSITION",MBVERTEX));
      hi_niit = row_multiIndex->get<Composite_mi_tag2>().upper_bound(boost::make_tuple("SPATIAL_POSITION",MBVERTEX));
      assert(distance(niit,hi_niit)==12);   
      for(int dd = 0;niit!=hi_niit;niit++,dd++) dofs_x[dd] = niit->get_FieldData(); 
      } catch (const char* msg) {
	SETERRQ(PETSC_COMM_SELF,1,msg);
      } 
    }
    break;
    default:
    SETERRQ(PETSC_COMM_SELF,1,"no implemented");
  }
  PetscFunctionReturn(0);
}
PetscErrorCode FEMethod_ComplexForLazy::GetTangent() {
  PetscFunctionBegin;
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
    double center[3]; 
    tetcircumcenter_tp(&coords[0],&coords[3],&coords[6],&coords[9],center,NULL,NULL,NULL); 
    double r = cblas_dnrm2(3,center,1);
    int g_dim = get_dim_gNTET();
    if(type_of_analysis&spatail_analysis) {
	assert(12 == RowGlob[0].size());
	Khh = ublas::zero_matrix<double>(12,12);
	assert(3*(unsigned int)NBVOLUME_H1(order_volume) == RowGlob[1+6+4].size());
	Kvolumeh.resize(RowGlob[1+6+4].size(),12);
	ee = 0;
	for(;ee<6;ee++) {
	  assert(3*(unsigned int)NBEDGE_H1(order_edges[ee]) == RowGlob[1+ee].size());
	  Kedgeh_data[ee].resize(RowGlob[1+ee].size(),12);
	  Kedgeh[ee] = &Kedgeh_data[ee].data()[0];
	  for(int eee = 0;eee<6;eee++) {
	    Khh_edgeedge_data(eee,ee).resize(RowGlob[1+eee].size(),RowGlob[1+ee].size());
	    Khh_edgeedge[eee][ee] = &Khh_edgeedge_data(eee,ee).data()[0];
	  }
	  for(int fff = 0;fff<4;fff++) {
	    Khh_faceedge_data(fff,ee).resize(RowGlob[1+6+fff].size(),RowGlob[1+ee].size());
	    Khh_faceedge[fff][ee] = &Khh_faceedge_data(fff,ee).data()[0];
	  }
	  Khedge_data[ee].resize(12,RowGlob[1+ee].size());
	  Khedge[ee] = &Khedge_data[ee].data()[0];
	  Khh_volumeedge_data[ee].resize(RowGlob[1+6+4].size(),RowGlob[1+ee].size());
	  Khh_volumeedge[ee] = & Khh_volumeedge_data[ee].data()[0];
	  //
	  Khh_edgevolume_data[ee].resize(RowGlob[1+ee].size(),RowGlob[1+6+4].size());
	  Khh_edgevolume[ee] = &Khh_edgevolume_data[ee].data()[0];
	}
	ff = 0;
	for(;ff<4;ff++) {
	  assert(3*(unsigned int)NBFACE_H1(order_faces[ff]) == RowGlob[1+6+ff].size());
	  Kfaceh_data[ff] = ublas::zero_matrix<double>(RowGlob[1+6+ff].size(),12);
	  Kfaceh[ff] = &Kfaceh_data[ff].data()[0];
	  Khface_data[ff].resize(12,RowGlob[1+6+ff].size());
	  Khface[ff] = &Khface_data[ff].data()[0];
	  Khh_volumeface_data[ff].resize(RowGlob[1+6+4].size(),RowGlob[1+6+ff].size());
	  Khh_volumeface[ff] = &Khface_data[ff].data()[0];
	  for(int fff = 0;fff<4;fff++) {
	    Khh_faceface_data(fff,ff).resize(RowGlob[1+6+fff].size(),RowGlob[1+6+ff].size());
	    Khh_faceface[fff][ff] = &Khh_faceface_data(fff,ff).data()[0];
	  }
	  for(int eee = 0;eee<6;eee++) {
	    Khh_edgeface_data(eee,ff).resize(RowGlob[1+eee].size(),RowGlob[1+6+ff].size());
	    Khh_edgeface[eee][ff] = &Khh_edgeface_data(eee,ff).data()[0];
	  }
	  //
	  Khh_facevolume_data[ff].resize(RowGlob[1+6+ff].size(),RowGlob[1+6+4].size());
	  Khh_facevolume[ff] = &Khh_facevolume_data[ff].data()[0];
	}
	Khvolume.resize(12,RowGlob[1+6+4].size());
	Khh_volumevolume.resize(RowGlob[1+6+4].size(),RowGlob[1+6+4].size());
    }
    unsigned int sub_analysis_type = (spatail_analysis|material_analysis)&type_of_analysis;
    switch(sub_analysis_type) {
      case spatail_analysis: {
	ierr = Tangent_hh_hierachical(&order_edges[0],&order_faces[0],order_volume,V,eps*r,lambda,mu,ptr_matctx,
	  &diffNTETinvJac[0],&diff_edgeNinvJac[0],&diff_faceNinvJac[0],&diff_volumeNinvJac[0], 
	  &coords[0],&dofs_x[0],&dofs_x_edge[0],&dofs_x_face[0],&dofs_x_volume[0], 
	  &Khh.data()[0],NULL,Kedgeh,Kfaceh,&Kvolumeh.data()[0],g_dim,g_TET_W); CHKERRQ(ierr);
	ierr = Tangent_hh_hierachical_edge(&order_edges[0],&order_faces[0],order_volume,V,eps*r,lambda,mu,ptr_matctx, 
	  &diffNTETinvJac[0],&diff_edgeNinvJac[0],&diff_faceNinvJac[0],&diff_volumeNinvJac[0], 
	  &coords[0],&dofs_x[0],&dofs_x_edge[0],&dofs_x_face[0],&dofs_x_volume[0], 
	  &Khedge[0],NULL,Khh_edgeedge,Khh_faceedge,Khh_volumeedge, 
	  g_dim,g_TET_W); CHKERRQ(ierr);
	ierr = Tangent_hh_hierachical_face(&order_edges[0],&order_faces[0],order_volume,V,eps*r,lambda,mu,ptr_matctx, 
	  &diffNTETinvJac[0],&diff_edgeNinvJac[0],&diff_faceNinvJac[0],&diff_volumeNinvJac[0], 
	  &coords[0],&dofs_x[0],&dofs_x_edge[0],&dofs_x_face[0],&dofs_x_volume[0], 
	  &Khface[0],NULL,Khh_edgeface,Khh_faceface,Khh_volumeface, 
	  g_dim,g_TET_W); CHKERRQ(ierr);
	ierr = Tangent_hh_hierachical_volume(&order_edges[0],&order_faces[0],order_volume,V,eps*r,lambda,mu,ptr_matctx, 
	  &diffNTETinvJac[0],&diff_edgeNinvJac[0],&diff_faceNinvJac[0],diff_volumeNinvJac, 
	  &coords[0],&dofs_x[0],&dofs_x_edge[0],&dofs_x_face[0],&dofs_x_volume[0], 
	  &Khvolume.data()[0],NULL,Khh_edgevolume,Khh_facevolume,&Khh_volumevolume.data()[0], 
	  g_dim,g_TET_W); CHKERRQ(ierr);
      }
      break;
      default:
	SETERRQ(PETSC_COMM_SELF,1,"no implemented");
    }
  }
  break;
  default:
    SETERRQ(PETSC_COMM_SELF,1,"no implemented");
  }
  PetscFunctionReturn(0);
}
PetscErrorCode FEMethod_ComplexForLazy::GetFint() {
  PetscFunctionBegin;
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
    if(type_of_analysis&spatail_analysis) {
      ee = 0;
      for(;ee<6;ee++) {
	Fint_h_edge_data[ee].resize(RowGlob[1+ee].size());
	Fint_h_edge[ee] = &Fint_h_edge_data[ee].data()[0];
      }
      ff = 0;
      for(;ff<4;ff++) {
	Fint_h_face_data[ff].resize(RowGlob[1+6+ff].size());
	Fint_h_face[ff] = &Fint_h_face_data[ff].data()[0];
      }
      Fint_h_volume.resize(RowGlob[1+6+4].size());
    }
    unsigned int sub_analysis_type = (spatail_analysis|material_analysis)&type_of_analysis;
    switch(sub_analysis_type) {
      case spatail_analysis: {
	ierr = Fint_Hh_hierarchical(&order_edges[0],&order_faces[0],order_volume,V,lambda,mu,ptr_matctx, 
	  &diffNTETinvJac[0],&diff_edgeNinvJac[0],&diff_faceNinvJac[0],&diff_volumeNinvJac[0], 
	  &coords[0],&dofs_x[0],NULL,NULL,
	  &dofs_x_edge[0],&dofs_x_face[0],&dofs_x_volume[0], 
	  NULL,&Fblock_x.data()[0],Fint_h_edge,Fint_h_face,&Fint_h_volume.data()[0],
	  NULL,NULL,NULL,NULL,NULL,
	  g_dim,g_TET_W); CHKERRQ(ierr);
      }
      break;
      default:
	SETERRQ(PETSC_COMM_SELF,1,"no implemented");
    }
  }
  break;
  default:
    SETERRQ(PETSC_COMM_SELF,1,"no implemented");
  }

  PetscFunctionReturn(0);
}
PetscErrorCode FEMethod_ComplexForLazy::GetFaceIndicesAndData(EntityHandle face) {
  PetscFunctionBegin;
  typedef FENumeredDofMoFEMEntity_multiIndex::index<Composite_mi_tag3>::type::iterator dofs_iterator;
  dofs_iterator fiit,hi_fiit;
  fiit = row_multiIndex->get<Composite_mi_tag3>().lower_bound(boost::make_tuple("SPATIAL_POSITION",face));
  if(fiit==row_multiIndex->get<Composite_mi_tag3>().end()) SETERRQ(PETSC_COMM_SELF,1,"no such ent");
  if(fiit->get_ent_type()!=MBTRI) SETERRQ(PETSC_COMM_SELF,1,"works only for facec");
  hi_fiit = row_multiIndex->get<Composite_mi_tag3>().upper_bound(boost::make_tuple("SPATIAL_POSITION",face));
  FaceIndices.resize(distance(fiit,hi_fiit));
  FaceData.resize(distance(fiit,hi_fiit));
  face_order = fiit->get_max_order();
  int dd = 0;
  if(NBFACE_H1(face_order)>0) {
    for(dofs_iterator fiiit = fiit;fiiit!=hi_fiit;fiiit++,dd++) {
      FaceIndices[dd] = fiiit->get_petsc_gloabl_dof_idx();
      FaceData[dd] = fiiit->get_FieldData();
    }
  }
  N_face.resize(g_TRI_dim*NBFACE_H1(face_order));
  diffN_face.resize(2*g_TRI_dim*NBFACE_H1(face_order));
  ierr = H1_FaceShapeFunctions_MBTRI(face_order,&g_NTRI[0],&diffNTRI[0],&N_face[0],&diffN_face[0],g_TRI_dim); CHKERRQ(ierr);
  NodeIndices.resize(12);
  NodeData.resize(12);
  const EntityHandle* conn_face; 
  int num_nodes; 
  rval = moab.get_connectivity(face,conn_face,num_nodes,true); CHKERR_PETSC(rval);
  int nn = 0;
  for(dd = 0;nn<4;nn++) {
    dofs_iterator niit,hi_niit;
    niit = row_multiIndex->get<Composite_mi_tag3>().lower_bound(boost::make_tuple("SPATIAL_POSITION",conn_face[nn]));
    hi_niit = row_multiIndex->get<Composite_mi_tag3>().upper_bound(boost::make_tuple("SPATIAL_POSITION",conn_face[nn]));
    for(;niit!=hi_niit;niit++,dd++) {
      NodeIndices[dd] = niit->get_petsc_gloabl_dof_idx();
      NodeData[dd] = niit->get_FieldData();
    }
  }
  EdgeIndices_data.resize(3);
  EdgeData_data.resize(3);
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
    eiit = row_multiIndex->get<Composite_mi_tag3>().lower_bound(boost::make_tuple("SPATIAL_POSITION",edge));
    hi_eiit = row_multiIndex->get<Composite_mi_tag3>().upper_bound(boost::make_tuple("SPATIAL_POSITION",edge));
    FaceEdgeOrder[ee] = eiit->get_max_order();
    if(NBEDGE_H1(FaceEdgeOrder[ee])>0) {
      EdgeIndices_data[ee].resize(distance(eiit,hi_eiit));
      EdgeData_data[ee].resize(distance(eiit,hi_eiit));
      for(dd = 0;eiit!=hi_eiit;eiit++,dd++) {
	EdgeIndices_data[ee][dd] = eiit->get_petsc_gloabl_dof_idx();
	EdgeData_data[ee][dd] = eiit->get_FieldData();
      }
      EdgeData[ee] = &(EdgeData_data[ee].data()[0]);
      N_edge_data[ee].resize(g_TRI_dim*NBEDGE_H1(FaceEdgeOrder[ee]));
      diffN_edge_data[ee].resize(2*g_TRI_dim*NBEDGE_H1(FaceEdgeOrder[ee]));
      N_edge[ee] = &(N_edge_data[ee][0]);
      diffN_edge[ee] = &(diffN_edge_data[ee][0]);
    }
  }
  ierr = H1_EdgeShapeFunctions_MBTRI(&FaceEdgeSense[0],&FaceEdgeOrder[0],&g_NTRI[0],diffNTRI,N_edge,diffN_edge,g_TRI_dim); CHKERRQ(ierr);
  GetFaceIndicesAndData_face = face;
  PetscFunctionReturn(0);
}
PetscErrorCode FEMethod_ComplexForLazy::GetFExt(EntityHandle face,double *t,double *t_edge[],double *t_face) {
  PetscFunctionBegin;
  if(GetFaceIndicesAndData_face!=face) SETERRQ(PETSC_COMM_SELF,1,"run GetFaceIndicesAndData(face) before call of this function");
  FExt.resize(9);
  FExt_edge_data.resize(3);
  int ee = 0;
  for(;ee<3;ee++) {
    FExt_edge_data[ee].resize(EdgeIndices_data[ee].size());
    FExt_edge[ee] = &FExt_edge_data[ee].data()[0];
  }
  FExt_face.resize(FaceIndices.size());
  ierr = Fext_h_hierarchical(
    face_order,&FaceEdgeOrder[0],//2
    &g_NTRI[0],&N_face[0],N_edge,&diffNTRI[0],&diffN_face[0],diffN_edge,//8
    t,t_edge,t_face,//11
    &NodeData.data()[0],EdgeData,&FaceData.data()[0],//14
    NULL,NULL,NULL,//17
    &FExt.data()[0],FExt_edge,&FExt_face.data()[0],//20
    NULL,NULL,NULL,//23
    g_TRI_dim,g_TRI_W); CHKERRQ(ierr);
  PetscFunctionReturn(0);
}
PetscErrorCode FEMethod_ComplexForLazy::GetTangentExt(EntityHandle face,double *t,double *t_edge[],double *t_face) {
  PetscFunctionBegin;
  if(GetFaceIndicesAndData_face!=face) SETERRQ(PETSC_COMM_SELF,1,"run GetFaceIndicesAndData(face) before call of this function");
  Kext_hh.resize(9,9);
  Kext_edgeh_data.resize(3);
  int ee = 0;
  for(;ee<3;ee++) {
    Kext_edgeh_data[ee].resize(EdgeIndices_data[ee].size(),9);
    Kext_edgeh[ee] = &Kext_edgeh_data[ee].data()[0];
  }
  Kext_faceh.resize(FaceIndices.size(),9);
  ierr = Kext_hh_hierarchical(eps,face_order,&FaceEdgeOrder[0],
    &gNTRI[0],&N_face[0],N_edge,&diffNTRI[0],&diffN_face[0],diffN_edge,
    t,t_edge,t_face,&NodeData.data()[0],EdgeData,&FaceData.data()[0],
    &Kext_hh.data()[0],Kext_edgeh,&Kext_faceh.data()[0],g_TRI_dim,g_TRI_W); CHKERRQ(ierr);
  Kext_hedge_data.resize(3);
  for(ee = 0;ee<3;ee++) {
    for(int eee = 0;eee<3;eee++) {
      Kext_hedge_data[ee].resize(EdgeIndices_data[ee].size(),EdgeIndices_data[eee].size());
    }

  }
  ierr = Kext_hh_hierarchical_edge(eps,face_order,&FaceEdgeOrder[0],
    &gNTRI[0],&N_face[0],N_edge,&diffNTRI[0],&diffN_face[0],diffN_edge,
    t,t_edge,t_face,&NodeData.data()[0],EdgeData,&FaceData.data()[0],
    Kext_hedge,Kext_edgeegde,Kext_faceedge,g_TRI_dim,g_TRI_W); CHKERRQ(ierr);


  /*PetscErrorCode Kext_hh_hierarchical_face(double eps,int order,int *order_edge,
  double *N,double *N_face,double *N_edge[],
  double *diffN,double *diffN_face,double *diffN_edge[],
  double *t,double *t_edge[],double *t_face,
  double *dofs_x,double *dofs_x_edge[],double *dofs_x_face,
  double *idofs_x_face,
  double *Kext_hface,double *Kext_egdeface[3],double *Kext_faceface,*/

  PetscFunctionReturn(0);
}


}

