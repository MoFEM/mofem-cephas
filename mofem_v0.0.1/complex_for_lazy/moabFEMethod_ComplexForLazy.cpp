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

namespace MoFEM {

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
      ee = 0;
      for(;ee<6;ee++) {
	FENumeredDofMoFEMEntity_multiIndex::index<Composite_mi_tag>::type::iterator eiit;
	eiit = row_multiIndex->get<Composite_mi_tag>().find(boost::make_tuple("SPATIAL_POSITION",MBEDGE,ff));
	if(eiit!=row_multiIndex->get<Composite_mi_tag>().end()) {
	  order_edges[ee] = eiit->get_max_order();
	} else {
	  order_edges[ee] = 0;
	}
      }
      ff = 0;
      for(;ff<4;ff++) {
	FENumeredDofMoFEMEntity_multiIndex::index<Composite_mi_tag>::type::iterator fiit;
	fiit = row_multiIndex->get<Composite_mi_tag>().find(boost::make_tuple("SPATIAL_POSITION",MBTRI,ff));
	if(fiit!=row_multiIndex->get<Composite_mi_tag>().end()) {
	  order_faces[ee] = fiit->get_max_order();
	} else {
	  order_faces[ee] = 0;
	}
      }
      FENumeredDofMoFEMEntity_multiIndex::index<Composite_mi_tag2>::type::iterator viit;
      viit = row_multiIndex->get<Composite_mi_tag2>().find(boost::make_tuple("SPATIAL_POSITION",MBTET));
      order_volume = viit->get_max_order();
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


    /*FENumeredDofMoFEMEntity_multiIndex::index<Composite_mi_tag>::type::iterator eiit,hi_eiit;
    eiit = row_multiIndex->get<Composite_mi_tag>().find(boost::make_tuple(field_name,MBEDGE,0));
    hi_eiit = row_multiIndex->get<Composite_mi_tag>().find(boost::make_tuple(field_name,MBEDGE,6));*/


    unsigned int sub_analysis_type = (spatail_analysis|material_analysis)&type_of_analysis;
    switch(sub_analysis_type) {
      case spatail_analysis:
	/*
	Tangent_hh_hierachical(order_edges,order_faces,order_volume,V,eps*r,lambda,mu,ptr_matctx, 
	  diffNinvJac,diff_edgeNinvJac,diff_faceNinvJac,diff_volumeNinvJac, 
	  dofs_X,dofs_x,dofs_x_edge,dofs_x_face,dofs_x_volume, 
	  Khh,NULL,Kedgeh,Kfaceh,Kvolumeh,g_dim,g_w); 
	Tangent_hh_hierachical_edge(order_edges,order_faces,order_volume,V,eps*r,lambda,mu,ptr_matctx, 
	  diffNinvJac,diff_edgeNinvJac,diff_faceNinvJac,diff_volumeNinvJac, 
	  dofs_X,dofs_x,dofs_x_edge,dofs_x_face,dofs_x_volume, 
	  Khedge,NULL,Khh_edgeedge,Khh_faceedge,Khh_volumeedge, 
	  g_dim,g_w); 
	Tangent_hh_hierachical_face(order_edges,order_faces,order_volume,V,eps*r,lambda,mu,ptr_matctx, 
	  diffNinvJac,diff_edgeNinvJac,diff_faceNinvJac,diff_volumeNinvJac, 
	  dofs_X,dofs_x,dofs_x_edge,dofs_x_face,dofs_x_volume, 
	  Khface,NULL,Khh_edgeface,Khh_faceface,Khh_volumeface, 
	  g_dim,g_w); 
	Tangent_hh_hierachical_volume(order_edges,order_faces,order_volume,V,eps*r,lambda,mu,ptr_matctx, 
	  diffNinvJac,diff_edgeNinvJac,diff_faceNinvJac,diff_volumeNinvJac, 
	  dofs_X,dofs_x,dofs_x_edge,dofs_x_face,dofs_x_volume, 
	  Khvolume,NULL,Khh_edgevolume,Khh_facevolume,Khh_volumevolume, 
	  g_dim,g_w);
	*/
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
  PetscFunctionReturn(0);
}
PetscErrorCode FEMethod_ComplexForLazy::GetFext() {
  PetscFunctionBegin;
  PetscFunctionReturn(0);
}

}

