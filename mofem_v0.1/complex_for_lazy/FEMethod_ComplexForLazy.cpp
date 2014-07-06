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

FEMethod_ComplexForLazy::FEMethod_ComplexForLazy(FieldInterface& _mField,analysis _type,double _lambda,double _mu,double _thermal_expansion,int _verbose):
    FEMethod_ComplexForLazy_Data(_mField,_verbose),
    type_of_analysis(_type),type_of_forces(conservative),
    lambda(_lambda),mu(_mu),thermal_expansion(_thermal_expansion),ptr_matctx(NULL),
    eps(1e-10),thermal_load_factor(1),
    spatial_field_name("SPATIAL_POSITION"),
    material_field_name("MESH_NODE_POSITIONS"),
    termal_field_name("TEMPERATURE") {
  order_x_edges.resize(6);
  order_x_faces.resize(4);
  order_X_edges.resize(6);
  order_X_faces.resize(4);
  diff_edgeNinvJac.resize(6);
  diff_faceNinvJac.resize(4);
  edgeN.resize(6);
  faceN.resize(4);
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
  iFint_h.resize(12);
  iFint_h_edge_data.resize(6);
  iFint_h_face_data.resize(4);
  //
  Fint_H.resize(12);
  iFint_H.resize(12);
  //
  g_NTET.resize(4*45);
  ShapeMBTET(&g_NTET[0],G_TET_X45,G_TET_Y45,G_TET_Z45,45);
  g_TET_W = G_TET_W45;

  propeties_from_BLOCKSET_MAT_ELASTICSET = false;
  for(_IT_CUBITMESHSETS_BY_BCDATA_TYPE_FOR_LOOP_(mField,BLOCKSET|MAT_ELASTICSET,it)) {
    propeties_from_BLOCKSET_MAT_ELASTICSET = true;
  }

}
FEMethod_ComplexForLazy::~FEMethod_ComplexForLazy() {
}
PetscErrorCode FEMethod_ComplexForLazy::GetMatParameters(double *_lambda,double *_mu,double *_thermal_expansion,void **ptr_matctx) {
  PetscFunctionBegin;

  *_lambda = lambda;
  *_mu = mu;
  *_thermal_expansion = thermal_expansion;

  if(propeties_from_BLOCKSET_MAT_ELASTICSET) {
    EntityHandle ent = fe_ptr->get_ent();
    for(_IT_CUBITMESHSETS_BY_BCDATA_TYPE_FOR_LOOP_(mField,BLOCKSET|MAT_ELASTICSET,it)) {

      Mat_Elastic mydata;
      ierr = it->get_attribute_data_structure(mydata); CHKERRQ(ierr);

      Range meshsets;
      rval = moab.get_entities_by_type(it->meshset,MBENTITYSET,meshsets,true); CHKERR_PETSC(rval);
      meshsets.insert(it->meshset);
      for(Range::iterator mit = meshsets.begin();mit != meshsets.end(); mit++) {
        if( moab.contains_entities(*mit,&ent,1) ) {
          *_lambda = LAMBDA(mydata.data.Young,mydata.data.Poisson);
          *_mu = MU(mydata.data.Young,mydata.data.Poisson);
          int material_type = (int)mydata.data.User1;
          if(it->get_Cubit_name().compare(0,29,"MAT_ELASTIC_EberleinHolzapfel") == 0) {
            Mat_Elastic_EberleinHolzapfel1 mydata;
            ierr = it->get_attribute_data_structure(mydata); CHKERRQ(ierr);
            set_PhysicalEquationNumber(eberleinholzapfel1);
            EberleinHolzapfel1_mat_parameters.eq_solid = neohookean;
            EberleinHolzapfel1_mat_parameters.k1 = mydata.data.k1;
            EberleinHolzapfel1_mat_parameters.k2 = mydata.data.k2;
            EberleinHolzapfel1_mat_parameters.fibre_vector_a1[0] = mydata.data.a0x;
            EberleinHolzapfel1_mat_parameters.fibre_vector_a1[1] = mydata.data.a0y;
            EberleinHolzapfel1_mat_parameters.fibre_vector_a1[2] = mydata.data.a0z;
            EberleinHolzapfel1_mat_parameters.fibre_vector_a2[0] = mydata.data.a1x;
            EberleinHolzapfel1_mat_parameters.fibre_vector_a2[1] = mydata.data.a1y;
            EberleinHolzapfel1_mat_parameters.fibre_vector_a2[2] = mydata.data.a1z;
            *ptr_matctx = &EberleinHolzapfel1_mat_parameters;
          } else {
            if(material_type>=10) {
              switch(material_type) {
                case 10: {
                           set_PhysicalEquationNumber(hooke);
                         }
                         break;
                case 11: {
                           set_PhysicalEquationNumber(stvenant_kirchhoff);
                         }
                         break;
                case 12: {
                           set_PhysicalEquationNumber(neohookean);
                         }
                         break;
                default: {
                           SETERRQ(PETSC_COMM_SELF,1,
                               "Materail not defined (Attribute 4):\n"
                               "\t10 = hooke\n"
                               "\t11 = stvenant_kirchhoff\n"
                               "\t12 = neohookean\n"
                               "\tname is MAT_ELASTIC_EberleinHolzapfel = eberleinholzapfel1\n");
                         }
              }
            }
          }
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
  vector<int>& order_edges,vector<int>& order_faces,int& order_volume,
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
	  if(diffH1edgeNinvJac[ee].size() < (unsigned int)NBEDGE_H1(order_edges[ee])) {
	    SETERRQ(PETSC_COMM_SELF,1,"not enugh shape functions calulated");
	  }
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
	  if(diffH1faceNinvJac[ff].size() < (unsigned int)NBFACE_H1(order_faces[ff])) {
	    SETERRQ(PETSC_COMM_SELF,1,"not enugh shape functions calulated");
	  }
	} else {
	  order_faces[ff] = 0;
	  dofs_face[ff] = 0;
	}
      }
      //data voolume
      FEDofMoFEMEntity_multiIndex::index<Composite_Name_And_Type_mi_tag>::type::iterator viit,hi_viit;
      viit = data_multiIndex->get<Composite_Name_And_Type_mi_tag>().lower_bound(boost::make_tuple(field_name,MBTET));
      hi_viit = data_multiIndex->get<Composite_Name_And_Type_mi_tag>().upper_bound(boost::make_tuple(field_name,MBTET));
      if(viit!=hi_viit) {
	order_volume = viit->get_max_order();
	dofs_volume.resize(distance(viit,hi_viit));
	assert(dofs_volume.size() == (unsigned int)3*NBVOLUME_H1(order_volume));
	for(;viit!=hi_viit;viit++) dofs_volume[viit->get_EntDofIdx()] = viit->get_FieldData();
	if(diffH1elemNinvJac.size() < (unsigned int)NBVOLUME_H1(order_volume)) {
	  SETERRQ(PETSC_COMM_SELF,1,"not enugh shape functions calulated");
	}
      } else {
	order_volume = 0;
      }
      //data nodes
      copy(coords.begin(),coords.end(),dofs_nodes.begin());
      FEDofMoFEMEntity_multiIndex::index<Composite_Name_Type_And_Side_Number_mi_tag>::type::iterator niit,hi_niit;
      niit = data_multiIndex->get<Composite_Name_Type_And_Side_Number_mi_tag>().lower_bound(boost::make_tuple(field_name,MBVERTEX,0));
      hi_niit = data_multiIndex->get<Composite_Name_Type_And_Side_Number_mi_tag>().upper_bound(boost::make_tuple(field_name,MBVERTEX,4));
      if(field_name == spatial_field_name) {
	if(distance(niit,hi_niit)!=12) {
	  EntityHandle ent = fe_ptr->get_ent();
	  const EntityHandle* conn;
	  int num_nodes;
	  rval = moab.get_connectivity(ent,conn,num_nodes,true); CHKERR_PETSC(rval);
	  PetscPrintf(PETSC_COMM_WORLD,"== sides ==\n");
	  SideNumber_multiIndex::nth_index<2>::type &sides =
	    const_cast<SideNumber_multiIndex::nth_index<2>::type&>(fe_ptr->get_side_number_table().get<2>());
	  SideNumber_multiIndex::nth_index<2>::type::iterator sit,hi_sit;
	  sit = sides.lower_bound(MBVERTEX);
	  hi_sit = sides.upper_bound(MBVERTEX);
	  for(;sit!=hi_sit;sit++) {
	    EntityHandle side_ent;
	    rval = moab.side_element(ent,0,sit->side_number,side_ent); CHKERR_PETSC(rval);
	    PetscPrintf(PETSC_COMM_WORLD,"side %u ent %lu (%lu) [%lu]\n",sit->side_number,sit->ent,conn[sit->side_number],side_ent);
	    DofMoFEMEntity_multiIndex::index<Composite_Name_And_Ent_mi_tag>::type::iterator dit,hi_dit;
	    dit = dofs_moabfield->get<Composite_Name_And_Ent_mi_tag>().lower_bound(boost::make_tuple(spatial_field_name,sit->ent));
	    hi_dit = dofs_moabfield->get<Composite_Name_And_Ent_mi_tag>().upper_bound(boost::make_tuple(spatial_field_name,sit->ent));
	    for(;dit!=hi_dit;dit++) {
	      ostringstream ss;
	      ss << *dit << endl;
	      PetscPrintf(PETSC_COMM_WORLD,"%s",ss.str().c_str());
	    }
	    MoFEMEntity_multiIndex::index<Composite_Name_And_Ent_mi_tag>::type::iterator eit;
	    eit = ents_moabfield->get<Composite_Name_And_Ent_mi_tag>().find(boost::make_tuple(spatial_field_name,sit->ent));
	    if(eit != ents_moabfield->get<Composite_Name_And_Ent_mi_tag>().end()) {
	      ostringstream ss;
	      ss << *eit << endl;
	      PetscPrintf(PETSC_COMM_WORLD,"%s",ss.str().c_str());
	    }
	    RefMoFEMEntity_multiIndex::index<MoABEnt_mi_tag>::type::iterator reit;
	    reit = refined_entities->get<MoABEnt_mi_tag>().find(sit->ent);
	    if(reit != refined_entities->get<MoABEnt_mi_tag>().end()) {
	      ostringstream ss;
	      ss << *reit << endl;
	      PetscPrintf(PETSC_COMM_WORLD,"%s",ss.str().c_str());
	    }
	  }
	  PetscPrintf(PETSC_COMM_WORLD,"== fe dofs ==\n");
	  FEDofMoFEMEntity_multiIndex::index<FieldName_mi_tag>::type::iterator iit,hi_iit;
	  iit = data_multiIndex->get<FieldName_mi_tag>().lower_bound(spatial_field_name);
	  hi_iit = data_multiIndex->get<FieldName_mi_tag>().upper_bound(spatial_field_name);
	  for(;iit!=hi_iit;iit++) {
	    ostringstream ss;
	    ss << *iit << endl;
	    PetscPrintf(PETSC_COMM_WORLD,"%s",ss.str().c_str());
	  }
	  PetscPrintf(PETSC_COMM_WORLD,"== fe ==\n");
	  {
	    ostringstream ss;
	    ss << *fe_ptr << endl;
	    PetscPrintf(PETSC_COMM_WORLD,"%s",ss.str().c_str());
	  }
	  SETERRQ1(PETSC_COMM_SELF,1,
	    "I can not find dofs on vertices, it should be 12 dofs (i.e. 4 nodes and 3 dofs for each node), but is %u",
	    distance(niit,hi_niit));
	}
      }
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
  ierr = GetData(
    order_x_edges,order_x_faces,order_x_volume,
    dofs_x_edge_data,dofs_x_edge,
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
PetscErrorCode FEMethod_ComplexForLazy::GetDofs_Termal_FromElementData() {
  PetscFunctionBegin;
  dofs_temp.resize(4);
  fill(dofs_temp.begin(),dofs_temp.end(),0);
  FEDofMoFEMEntity_multiIndex::index<Composite_Name_Type_And_Side_Number_mi_tag>::type::iterator niit,hi_niit;
  niit = data_multiIndex->get<Composite_Name_Type_And_Side_Number_mi_tag>().lower_bound(boost::make_tuple(termal_field_name,MBVERTEX,0));
  hi_niit = data_multiIndex->get<Composite_Name_Type_And_Side_Number_mi_tag>().upper_bound(boost::make_tuple(termal_field_name,MBVERTEX,4));
  for(;niit!=hi_niit;niit++) {
    dofs_temp[niit->side_number_ptr->side_number] = niit->get_FieldData();
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
      edgeN[ee] = &*H1edgeN[ee].begin();
    }
    int ff = 0;
    for(;ff<4;ff++) {
      diff_faceNinvJac[ff] = &*diffH1faceNinvJac[ff].begin();
      faceN[ff] = &*H1faceN[ff].begin();
    }
    diff_volumeNinvJac = &*diffH1elemNinvJac.begin();
    volumeN = &*H1elemN.begin();
    double center[3];
    tetcircumcenter_tp(&coords[0],&coords[3],&coords[6],&coords[9],center,NULL,NULL,NULL);
    cblas_daxpy(3,-1,&coords[0],1,center,1);
    double r = cblas_dnrm2(3,center,1);
    int g_dim = get_dim_gNTET();
    if(type_of_analysis&spatail_analysis) {
	if(RowGlobSpatial.empty()) {
	  SETERRQ(PETSC_COMM_SELF,1,"RowGlobSpatial is not set");
	}
	if(12 != RowGlobSpatial[0].size()) {
	  SETERRQ(PETSC_COMM_SELF,1,"12 != RowGlobSpatial[0].size()");
	}
	KHh.resize(12,12);
	Khh.resize(12,12);
	ee = 0;
	for(;ee<6;ee++) {
	  assert(3*(unsigned int)NBEDGE_H1(order_x_edges[ee]) == RowGlobSpatial[1+ee].size());
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
	  assert(3*(unsigned int)NBFACE_H1(order_x_faces[ff]) == RowGlobSpatial[1+6+ff].size());
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
	assert(3*(unsigned int)NBVOLUME_H1(order_x_volume) == RowGlobSpatial[i_volume].size());
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
	  order_x_edges[ee] = 0;
	} else {
	  order_x_edges[ee] = eiit->get_max_order();
	}
	if(3*(unsigned int)NBEDGE_H1(order_x_edges[ee]) != dofs_x_edge_data[ee].size()) {
	  SETERRQ2(PETSC_COMM_SELF,1,"data inconsistency 3*NBEDGE_H1(order_x_edges[ee]) = %d != dofs_x_edge_data[ee].size() = %d ",
	    3*NBEDGE_H1(order_x_edges[ee]),dofs_x_edge_data[ee].size());
	}
	KedgeH_data[ee].resize(dofs_x_edge_data[ee].size(),12);
	KedgeH[ee] = &*KedgeH_data[ee].data().begin();
      }
      ff = 0;
      for(;ff<4;ff++) {
	FEDofMoFEMEntity_multiIndex::index<Composite_Name_Type_And_Side_Number_mi_tag>::type::iterator fiit;
	fiit = data_multiIndex->get<Composite_Name_Type_And_Side_Number_mi_tag>().find(boost::make_tuple(spatial_field_name,MBTRI,ff));
	if(fiit==data_multiIndex->get<Composite_Name_Type_And_Side_Number_mi_tag>().end()) {
	  order_x_faces[ff] = 0;
	} else {
	  order_x_faces[ff] = fiit->get_max_order();
	}
	assert(3*(unsigned int)NBFACE_H1(order_x_faces[ff]) == dofs_x_face_data[ff].size());
	KfaceH_data[ff].resize(dofs_x_face_data[ff].size(),12);
	KfaceH[ff] = &*KfaceH_data[ff].data().begin();
      }
      FEDofMoFEMEntity_multiIndex::index<Composite_Name_And_Type_mi_tag>::type::iterator viit;
      viit = data_multiIndex->get<Composite_Name_And_Type_mi_tag>().find(boost::make_tuple(spatial_field_name,MBTET));
      if(viit==data_multiIndex->get<Composite_Name_And_Type_mi_tag>().end()) {
	order_x_volume = 0;
      } else {
	order_x_volume = viit->get_max_order();
      }
      KvolumeH.resize(dofs_x_volume.size(),12);
    }
    ierr = GetData(
      order_X_edges,order_X_faces,order_X_volume,
      dofs_X_edge_data,dofs_X_edge,
      dofs_X_face_data,dofs_X_face,
      dofs_X_volume,dofs_X,
      material_field_name); CHKERRQ(ierr);
    ierr = GetDofs_Termal_FromElementData(); CHKERRQ(ierr);
    unsigned int sub_analysis_type[3] = {
      spatail_analysis, material_analysis, mesh_quality_analysis };
    double _lambda,_mu,_thermal_expansion;
    if( (spatail_analysis|material_analysis)&type_of_analysis ) {
      ierr = GetMatParameters(&_lambda,&_mu,&_thermal_expansion,&ptr_matctx); CHKERRQ(ierr);
    }
    int _order_T_volume = 0;
    for(int ss = 0;ss<3;ss++) {
      switch(sub_analysis_type[ss]&type_of_analysis) {
	case spatail_analysis: {
	  ierr = Tangent_hh_hierachical(
	    &maxOrderEdgeH1[0],&maxOrderFaceH1[0],maxOrderElemH1,
	    &order_X_edges[0],&order_X_faces[0],order_X_volume,
	    &order_x_edges[0],&order_x_faces[0],order_x_volume,
	    V,eps*r,_lambda,_mu,ptr_matctx,
	    &diffNTETinvJac[0],&diff_edgeNinvJac[0],&diff_faceNinvJac[0],&diff_volumeNinvJac[0],
	    &dofs_X[0],&dofs_X_edge[0],&dofs_X_face[0],&*dofs_X_volume.data().begin(),
	    &dofs_x[0],&dofs_x_edge[0],&dofs_x_face[0],&*dofs_x_volume.data().begin(),
	    //temperature
	    _thermal_expansion,thermal_load_factor,
	    &g_NTET[0],&edgeN[0],&faceN[0],volumeN, //shape functions
	    NULL,NULL,_order_T_volume, &dofs_temp.data()[0],NULL,NULL,NULL,
	    //
	    &*Khh.data().begin(),&*KHh.data().begin(),Kedgeh,Kfaceh,&*Kvolumeh.data().begin(),
	    g_dim,g_TET_W); CHKERRQ(ierr);
	  ierr = Tangent_hh_hierachical_edge(
	    &maxOrderEdgeH1[0],&maxOrderFaceH1[0],maxOrderElemH1,
	    &order_X_edges[0],&order_X_faces[0],order_X_volume,
	    &order_x_edges[0],&order_x_faces[0],order_x_volume,
	    V,eps*r,_lambda,_mu,ptr_matctx,
	    &diffNTETinvJac[0],&diff_edgeNinvJac[0],&diff_faceNinvJac[0],&diff_volumeNinvJac[0],
	    &dofs_X[0],&dofs_X_edge[0],&dofs_X_face[0],&*dofs_X_volume.data().begin(),
	    &dofs_x[0],&dofs_x_edge[0],&dofs_x_face[0],&*dofs_x_volume.data().begin(),
	    //temperature
	    _thermal_expansion,thermal_load_factor,
	    &g_NTET[0],&edgeN[0],&faceN[0],volumeN, //shape functions
	    NULL,NULL,_order_T_volume, &dofs_temp.data()[0],NULL,NULL,NULL,
	    //
	    &Khedge[0],&KHedge[0],Khh_edgeedge,Khh_faceedge,Khh_volumeedge,
	    g_dim,g_TET_W); CHKERRQ(ierr);
	  ierr = Tangent_hh_hierachical_face(
	    &maxOrderEdgeH1[0],&maxOrderFaceH1[0],maxOrderElemH1,
	    &order_X_edges[0],&order_X_faces[0],order_X_volume,
	    &order_x_edges[0],&order_x_faces[0],order_x_volume,
	    V,eps*r,_lambda,_mu,ptr_matctx,
	    &diffNTETinvJac[0],&diff_edgeNinvJac[0],&diff_faceNinvJac[0],&diff_volumeNinvJac[0],
	    &dofs_X[0],&dofs_X_edge[0],&dofs_X_face[0],&*dofs_X_volume.data().begin(),
	    &dofs_x[0],&dofs_x_edge[0],&dofs_x_face[0],&*dofs_x_volume.data().begin(),
	    //temperature
	    _thermal_expansion,thermal_load_factor,
	    &g_NTET[0],&edgeN[0],&faceN[0],volumeN, //shape functions
	    NULL,NULL,_order_T_volume, &dofs_temp.data()[0],NULL,NULL,NULL,
	    //
	    &Khface[0],&KHface[0],Khh_edgeface,Khh_faceface,Khh_volumeface,
	    g_dim,g_TET_W); CHKERRQ(ierr);
	  ierr = Tangent_hh_hierachical_volume(
	    &maxOrderEdgeH1[0],&maxOrderFaceH1[0],maxOrderElemH1,
	    &order_X_edges[0],&order_X_faces[0],order_X_volume,
	    &order_x_edges[0],&order_x_faces[0],order_x_volume,
	    V,eps*r,_lambda,_mu,ptr_matctx,
	    &diffNTETinvJac[0],&diff_edgeNinvJac[0],&diff_faceNinvJac[0],diff_volumeNinvJac,
	    &dofs_X[0],&dofs_X_edge[0],&dofs_X_face[0],&*dofs_X_volume.data().begin(),
	    &dofs_x[0],&dofs_x_edge[0],&dofs_x_face[0],&*dofs_x_volume.data().begin(),
	    //temperature
	    _thermal_expansion,thermal_load_factor,
	    &g_NTET[0],&edgeN[0],&faceN[0],volumeN, //shape functions
	    NULL,NULL,_order_T_volume, &dofs_temp.data()[0],NULL,NULL,NULL,
	    //
	    &*Khvolume.data().begin(),&*KHvolume.data().begin(),Khh_edgevolume,Khh_facevolume,&*Khh_volumevolume.data().begin(),
	    g_dim,g_TET_W); CHKERRQ(ierr);
	}
	break;
	case material_analysis: {
	  ierr = Tangent_HH_hierachical(
	    &maxOrderEdgeH1[0],&maxOrderFaceH1[0],maxOrderElemH1,
	    &order_X_edges[0],&order_X_faces[0],order_X_volume,
	    &order_x_edges[0],&order_x_faces[0],order_x_volume,
	    V,eps*r,_lambda,_mu,ptr_matctx,
	    &diffNTETinvJac[0],&diff_edgeNinvJac[0],&diff_faceNinvJac[0],&diff_volumeNinvJac[0],
	    &dofs_X[0],&dofs_X_edge[0],&dofs_X_face[0],&*dofs_X_volume.data().begin(),
	    &dofs_x[0],&dofs_x_edge[0],&dofs_x_face[0],&*dofs_x_volume.data().begin(),
	    //temperature
	    _thermal_expansion,thermal_load_factor,
	    &g_NTET[0],&edgeN[0],&faceN[0],volumeN, //shape functions
	    NULL,NULL,_order_T_volume, &dofs_temp.data()[0],NULL,NULL,NULL,
	    //
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
PetscErrorCode FEMethod_ComplexForLazy::GetFint() {
  PetscFunctionBegin;
  try {
  switch (fe_ent_ptr->get_ent_type()) {
    case MBTET: {
      int ee = 0;
      for(;ee<6;ee++) {
  	diff_edgeNinvJac[ee] = &(diffH1edgeNinvJac[ee])[0];
	edgeN[ee] = &*H1edgeN[ee].begin();
      }
      int ff = 0;
      for(;ff<4;ff++) {
        diff_faceNinvJac[ff] = &(diffH1faceNinvJac[ff])[0];
	faceN[ff] = &*H1faceN[ff].begin();
      }
      diff_volumeNinvJac = &diffH1elemNinvJac[0];
      volumeN = &*H1elemN.begin();
      int g_dim = get_dim_gNTET();
      if(!dofs_x_edge_data.empty()) {
	if(dofs_x_edge_data.size()!=6) SETERRQ(PETSC_COMM_SELF,1,"data vectors are not set");
	ee = 0;
	for(;ee<6;ee++) {
	  if(dofs_x_edge_data[ee].size()!=0) {
	    Fint_h_edge_data[ee].resize(dofs_x_edge_data[ee].size());
	    Fint_h_edge[ee] = &Fint_h_edge_data[ee].data()[0];
	    iFint_h_edge_data[ee].resize(dofs_x_edge_data[ee].size());
	    iFint_h_edge[ee] = &iFint_h_edge_data[ee].data()[0];
	  } else {
	    Fint_h_edge[ee] = NULL;
	    iFint_h_edge[ee] = NULL;
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
	    iFint_h_face_data[ff].resize(dofs_x_face_data[ff].size());
	    iFint_h_face[ff] = &iFint_h_face_data[ff].data()[0];
	  } else {
	    Fint_h_face[ff] = NULL;
	    iFint_h_face[ff] = NULL;
	  }
	}
      }
      if(dofs_x_volume.size()!=0) {
	  Fint_h_volume.resize(dofs_x_volume.size());
	  iFint_h_volume.resize(dofs_x_volume.size());
      }
      ierr = GetData(
	order_X_edges,order_X_faces,order_X_volume,
	dofs_X_edge_data,dofs_X_edge,
	dofs_X_face_data,dofs_X_face,
	dofs_X_volume,dofs_X,
	material_field_name); CHKERRQ(ierr);
      ierr = GetDofs_Termal_FromElementData(); CHKERRQ(ierr);
      unsigned int sub_analysis_type[5] = {
	spatail_analysis, material_analysis, mesh_quality_analysis, scaled_themp_direvative_spatial, scaled_themp_direvative_material
      };
      double _lambda,_mu,_thermal_expansion;
      if( (spatail_analysis|material_analysis|scaled_themp_direvative_spatial|scaled_themp_direvative_material)&type_of_analysis ) {
	ierr = GetMatParameters(&_lambda,&_mu,&_thermal_expansion,&ptr_matctx); CHKERRQ(ierr);
      }
      int _order_T_volume = 0;
      for(int ss = 0;ss<5;ss++) {
	switch(sub_analysis_type[ss]&type_of_analysis) {
	  case spatail_analysis: {
	    ierr = Fint_Hh_hierarchical(
	      &maxOrderEdgeH1[0],&maxOrderFaceH1[0],maxOrderElemH1,
	      &order_X_edges[0],&order_X_faces[0],order_X_volume,
	      &order_x_edges[0],&order_x_faces[0],order_x_volume,
	      V,_lambda,_mu,ptr_matctx,
	      &diffNTETinvJac[0],&diff_edgeNinvJac[0],&diff_faceNinvJac[0],&diff_volumeNinvJac[0],
	      &dofs_X.data()[0],&dofs_X_edge[0],&dofs_X_face[0],&*dofs_X_volume.data().begin(),
	      &*dofs_x.data().begin(),&dofs_x_edge[0],&dofs_x_face[0],&*dofs_x_volume.data().begin(),
	      NULL,NULL,
	      //temperature
	      _thermal_expansion,thermal_load_factor,eps,
	      &g_NTET[0],&edgeN[0],&faceN[0],volumeN, //shape functions
	      NULL,NULL,_order_T_volume, &dofs_temp.data()[0],NULL,NULL,NULL,
	      //
	      NULL,NULL,NULL,NULL,
	      &*Fint_h.data().begin(),Fint_h_edge,Fint_h_face,&*Fint_h_volume.data().begin(),
	      NULL,NULL,NULL,NULL,
	      NULL,NULL,NULL,NULL,
	      g_dim,g_TET_W); CHKERRQ(ierr);
	  }
	  break;
	  case material_analysis: {
	    ierr = Fint_Hh_hierarchical(
	      &maxOrderEdgeH1[0],&maxOrderFaceH1[0],maxOrderElemH1,
	      &order_X_edges[0],&order_X_faces[0],order_X_volume,
	      &order_x_edges[0],&order_x_faces[0],order_x_volume,
	      V,_lambda,_mu,ptr_matctx,
	      &diffNTETinvJac[0],&diff_edgeNinvJac[0],&diff_faceNinvJac[0],&diff_volumeNinvJac[0],
	      &dofs_X.data()[0],&dofs_X_edge[0],&dofs_X_face[0],&*dofs_X_volume.data().begin(),
	      &*dofs_x.data().begin(),&dofs_x_edge[0],&dofs_x_face[0],&*dofs_x_volume.data().begin(),
	      NULL,NULL,
	      //temperature
	      _thermal_expansion,thermal_load_factor,0,
	      &g_NTET[0],&edgeN[0],&faceN[0],volumeN, //shape functions
	      NULL,NULL,_order_T_volume, &dofs_temp.data()[0],NULL,NULL,NULL,
	      //
	      &*Fint_H.data().begin(),NULL,NULL,NULL,
	      NULL,NULL,NULL,NULL,
	      NULL,NULL,NULL,NULL,
	      NULL,NULL,NULL,NULL,
	      g_dim,g_TET_W); CHKERRQ(ierr);
	  }
	  break;
	  case mesh_quality_analysis: {
	    if(!flg_alpha2) SETERRQ(PETSC_COMM_SELF,1,"-my_alpha2 is not set");
	    if(!flg_gamma) SETERRQ(PETSC_COMM_SELF,1,"-my_gamma is not set");
	    if(*max_element(order_X_edges.begin(),order_X_edges.end())!=0) {
	      SETERRQ(PETSC_COMM_SELF,1,"not working for higher order tets");
	    }
	    if(*max_element(order_X_faces.begin(),order_X_faces.end())!=0) {
	      SETERRQ(PETSC_COMM_SELF,1,"not working for higher order tets");
	    }
	    if(order_X_volume!=0) {
	      SETERRQ(PETSC_COMM_SELF,1,"not working for higher order tets");
	    }
	    double coords_edges[2*3*6];
	    double alpha2_array[4] = { alpha2,alpha2,alpha2,alpha2 };
	    ierr = get_edges_from_elem_coords(&coords[0],coords_edges); CHKERRQ(ierr);
	    ierr = quality_volume_length_F(V,alpha2_array,gamma,&diffNTETinvJac[0],coords_edges,
	      &dofs_X.data()[0],NULL,NULL,NULL,quality0,quality,b,
	      &*Fint_H.data().begin(),NULL); CHKERRQ(ierr);
	  }
	  break;
	  case scaled_themp_direvative_spatial: {
	    ierr = Fint_Hh_hierarchical(
	      &maxOrderEdgeH1[0],&maxOrderFaceH1[0],maxOrderElemH1,
	      &order_X_edges[0],&order_X_faces[0],order_X_volume,
	      &order_x_edges[0],&order_x_faces[0],order_x_volume,
	      V,_lambda,_mu,ptr_matctx,
	      &diffNTETinvJac[0],&diff_edgeNinvJac[0],&diff_faceNinvJac[0],&diff_volumeNinvJac[0],
	      &dofs_X.data()[0],&dofs_X_edge[0],&dofs_X_face[0],&*dofs_X_volume.data().begin(),
	      &*dofs_x.data().begin(),&dofs_x_edge[0],&dofs_x_face[0],&*dofs_x_volume.data().begin(),
	      NULL,NULL,
	      //temperature
	      _thermal_expansion,thermal_load_factor,eps,
	      &g_NTET[0],&edgeN[0],&faceN[0],volumeN, //shape functions
	      NULL,NULL,_order_T_volume,
	      &dofs_temp.data()[0],NULL,NULL,NULL,
	      //
	      NULL,NULL,NULL,NULL,
	      NULL,NULL,NULL,NULL,
	      NULL,NULL,NULL,NULL,
	      &*iFint_h.data().begin(),iFint_h_edge,iFint_h_face,&*iFint_h_volume.data().begin(),
	      g_dim,g_TET_W); CHKERRQ(ierr);
	  }
	  break;
	  case scaled_themp_direvative_material: {
	    ierr = Fint_Hh_hierarchical(
	      &maxOrderEdgeH1[0],&maxOrderFaceH1[0],maxOrderElemH1,
	      &order_X_edges[0],&order_X_faces[0],order_X_volume,
	      &order_x_edges[0],&order_x_faces[0],order_x_volume,
	      V,_lambda,_mu,ptr_matctx,
	      &diffNTETinvJac[0],&diff_edgeNinvJac[0],&diff_faceNinvJac[0],&diff_volumeNinvJac[0],
	      &dofs_X.data()[0],&dofs_X_edge[0],&dofs_X_face[0],&*dofs_X_volume.data().begin(),
	      &*dofs_x.data().begin(),&dofs_x_edge[0],&dofs_x_face[0],&*dofs_x_volume.data().begin(),
	      NULL,NULL,
	      //temperature
	      _thermal_expansion,thermal_load_factor,eps,
	      &g_NTET[0],&edgeN[0],&faceN[0],volumeN, //shape functions
	      NULL,NULL,_order_T_volume,
	      &dofs_temp.data()[0],NULL,NULL,NULL,
	      //
	      NULL,NULL,NULL,NULL,
	      NULL,NULL,NULL,NULL,
	      &*iFint_H.data().begin(),NULL,NULL,NULL,
	      NULL,NULL,NULL,NULL,
	      g_dim,g_TET_W); CHKERRQ(ierr);
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


}

