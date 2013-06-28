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

#include<moabFEMethod_UpLevelStudent.hpp>
#include<FEM.h>

namespace MoFEM {

FEMethod_UpLevelStudent::FEMethod_UpLevelStudent(Interface& _moab,int _verbose): FEMethod_LowLevelStudent(_moab,_verbose) {
  double def_V = 0;
  rval = moab.tag_get_handle("Volume",1,MB_TYPE_DOUBLE,th_volume,MB_TAG_CREAT|MB_TAG_SPARSE,&def_V); CHKERR(rval);
}
FEMethod_UpLevelStudent::~FEMethod_UpLevelStudent() {}
PetscErrorCode FEMethod_UpLevelStudent::OpStudentStart_TET(vector<double>& _gNTET_) {
  PetscFunctionBegin;
    fe_ent_ptr = fe_ptr->fe_ptr;
    ierr = InitDataStructures(); CHKERRQ(ierr);
    ierr = GlobIndices(); CHKERRQ(ierr);
    ierr = DataOp(); CHKERRQ(ierr);
    ierr = ShapeFunctions_TET(_gNTET_); CHKERRQ(ierr);
    ierr = Data_at_GaussPoints(); CHKERRQ(ierr);
    ierr = DiffData_at_GaussPoints(); CHKERRQ(ierr);
    ierr = GetRowNMatrix_at_GaussPoint(); CHKERRQ(ierr);
    ierr = GetColNMatrix_at_GaussPoint(); CHKERRQ(ierr);
  try {
    ierr = GetRowDiffNMatrix_at_GaussPoint(); CHKERRQ(ierr);
    ierr = GetColDiffNMatrix_at_GaussPoint(); CHKERRQ(ierr);
  } catch (const std::exception& ex) {
    ostringstream ss;
    ss << "thorw in GetRowDiffNMatrix_at_GaussPoint(): " << ex.what() << endl;
    SETERRQ(PETSC_COMM_SELF,1,ss.str().c_str());
  }

  EntityHandle fe_handle = fe_ptr->get_ent();
  V = Shape_intVolumeMBTET(diffNTET,&coords[0]); 
  if( V <= 0 ) SETERRQ1(PETSC_COMM_SELF,1,"V < 0 for EntityHandle = %lu\n",fe_handle);
  rval = moab.tag_set_data(th_volume,&fe_handle,1,&V); CHKERR_PETSC(rval);

  const int g_dim = get_dim_gNTET();
  coords_at_Gauss_nodes.resize(g_dim);
  for(int gg = 0;gg<g_dim;gg++) {
    coords_at_Gauss_nodes[gg].resize(3);
    for(int dd = 0;dd<3;dd++) {
      (coords_at_Gauss_nodes[gg])[0] = cblas_ddot(4,&coords[0],3,&get_gNTET()[gg*4],1);
      (coords_at_Gauss_nodes[gg])[1] = cblas_ddot(4,&coords[1],3,&get_gNTET()[gg*4],1);
      (coords_at_Gauss_nodes[gg])[2] = cblas_ddot(4,&coords[2],3,&get_gNTET()[gg*4],1);
    }
  }

  PetscFunctionReturn(0);
}
PetscErrorCode FEMethod_UpLevelStudent::OpStudentStart_PRISM(vector<double>& _gNTRI_) {
  PetscFunctionBegin;
  fe_ent_ptr = fe_ptr->fe_ptr;
  ierr = InitDataStructures(); CHKERRQ(ierr);
  ierr = GlobIndices(); CHKERRQ(ierr);
  ierr = DataOp(); CHKERRQ(ierr);

  ierr = ShapeFunctions_PRISM(_gNTRI_); CHKERRQ(ierr);
  ierr = Data_at_GaussPoints(); CHKERRQ(ierr);
  ierr = GetRowNMatrix_at_GaussPoint(); CHKERRQ(ierr);
  ierr = GetColNMatrix_at_GaussPoint(); CHKERRQ(ierr);

  SideNumber_multiIndex& side_table = const_cast<SideNumber_multiIndex&>(fe_ent_ptr->get_side_number_table());
  SideNumber_multiIndex::nth_index<1>::type::iterator siit3 = side_table.get<1>().find(boost::make_tuple(MBTRI,3));
  SideNumber_multiIndex::nth_index<1>::type::iterator siit4 = side_table.get<1>().find(boost::make_tuple(MBTRI,4));
	
  ierr = ShapeDiffMBTRI(diffNTRI); CHKERRQ(ierr);

  int num_nodes;
  const EntityHandle *conn_prism;
  rval = moab.get_connectivity(fe_ent_ptr->get_ent(),conn_prism,num_nodes,true); CHKERR_PETSC(rval);
  assert(num_nodes==6);
  for(int nn = 0;nn<3;nn++) {
    conn_face3[nn] = conn_prism[nn];
    conn_face4[nn] = conn_prism[nn+3];
    //cerr << conn_face3[nn] << " ::: " << conn_face4[nn] << endl;
  }

  //face3
  rval = moab.get_coords(conn_face3,3,coords_face3); CHKERR_PETSC(rval);
  ierr = ShapeFaceNormalMBTRI(diffNTRI,coords_face3,normal3); CHKERRQ(ierr);
  area3 = cblas_dnrm2(3,normal3,1);
  //face4
  rval = moab.get_coords(conn_face4,3,coords_face4); CHKERR_PETSC(rval);
  ierr = ShapeFaceNormalMBTRI(diffNTRI,coords_face4,normal4); CHKERRQ(ierr);
  area4 = cblas_dnrm2(3,normal4,1);

  /*copy(&normal3[0],&normal3[3],ostream_iterator<double>(cerr," "));
  cerr << " -->  ";
  copy(&normal4[0],&normal4[3],ostream_iterator<double>(cerr," "));
  cerr << " : " << cblas_ddot(3,normal3,1,normal4,1)/(area3*area4) <<  endl;*/

  PetscFunctionReturn(0);
}
PetscErrorCode FEMethod_UpLevelStudent::OpStudentEnd() {
  PetscFunctionBegin;
  PetscFunctionReturn(0);
}
PetscErrorCode FEMethod_UpLevelStudent::GetRowIndices(const string &field_name,vector<DofIdx> &RowGlobDofs) {
  PetscFunctionBegin;
  MoFEMField_multiIndex::index<FieldName_mi_tag>::type::iterator fiit = moabfields->get<FieldName_mi_tag>().find(field_name);
  if(fiit==moabfields->get<FieldName_mi_tag>().end()) SETERRQ(PETSC_COMM_SELF,1,"no such field");
  GlobIndices_Type::iterator miit = row_nodesGlobIndices.find(fiit->get_MoFEMField_ptr());
  if(miit == row_nodesGlobIndices.end()) SETERRQ(PETSC_COMM_SELF,1,"no such field in FE");
  RowGlobDofs = miit->second;
  PetscFunctionReturn(0);
}
PetscErrorCode FEMethod_UpLevelStudent::GetRowIndices(const string &field_name,EntityType type,vector<DofIdx> &RowGlobDofs,int side_number) {
  PetscFunctionBegin;
  switch(type) {
    case MBEDGE:
    case MBTRI: {
      FENumeredDofMoFEMEntity_multiIndex::index<Composite_mi_tag>::type::iterator eiit;
      eiit = row_multiIndex->get<Composite_mi_tag>().find(boost::make_tuple(field_name,type,side_number));
      if(eiit == row_multiIndex->get<Composite_mi_tag>().end()) PetscFunctionReturn(0); //SETERRQ(PETSC_COMM_SELF,1,"no such ent");
      switch(type) {
	case MBEDGE: {
	  GlobIndices_EntType::iterator miit = row_edgesGlobIndices.find(eiit->get_MoFEMEntity_ptr());
	  if(miit == row_edgesGlobIndices.end()) SETERRQ(PETSC_COMM_SELF,1,"no such ent in FE");
	  RowGlobDofs = miit->second;
	} break;
	case MBTRI: {
	  GlobIndices_EntType::iterator miit = row_facesGlobIndices.find(eiit->get_MoFEMEntity_ptr());
	  if(miit == row_facesGlobIndices.end()) SETERRQ(PETSC_COMM_SELF,1,"no such ent in FE");
	  RowGlobDofs = miit->second;
	} break;
	default:
	SETERRQ(PETSC_COMM_SELF,1,"no implemented");
      }
    } break;
    case MBTET:
    case MBPRISM: {
      FENumeredDofMoFEMEntity_multiIndex::index<MoABEnt_mi_tag>::type::iterator eiit;
      eiit = row_multiIndex->get<MoABEnt_mi_tag>().find(fe_ent_ptr->get_ent());
      if(eiit == row_multiIndex->get<MoABEnt_mi_tag>().end()) PetscFunctionReturn(0);//SETERRQ(PETSC_COMM_SELF,1,"no such ent");
      GlobIndices_EntType::iterator miit = row_elemGlobIndices.find(eiit->get_MoFEMEntity_ptr());
      if(miit == row_elemGlobIndices.end()) SETERRQ(PETSC_COMM_SELF,1,"no such ent in FE");
      RowGlobDofs = miit->second;
    } break;
    default:
      SETERRQ(PETSC_COMM_SELF,1,"no implemented");
  }
  PetscFunctionReturn(0);
}
PetscErrorCode FEMethod_UpLevelStudent::GetColIndices(const string &field_name,vector<DofIdx> &ColGlobDofs) {
  PetscFunctionBegin;
  MoFEMField_multiIndex::index<FieldName_mi_tag>::type::iterator fiit = moabfields->get<FieldName_mi_tag>().find(field_name);
  if(fiit==moabfields->get<FieldName_mi_tag>().end()) SETERRQ(PETSC_COMM_SELF,1,"no such field");
  GlobIndices_Type::iterator miit = col_nodesGlobIndices.find(fiit->get_MoFEMField_ptr());
  if(miit == col_nodesGlobIndices.end()) SETERRQ(PETSC_COMM_SELF,1,"no such field in FE");
  ColGlobDofs = miit->second;
  PetscFunctionReturn(0);
}
PetscErrorCode FEMethod_UpLevelStudent::GetColIndices(const string &field_name,EntityType type,vector<DofIdx> &ColGlobDofs,int side_number) {
  PetscFunctionBegin;
  switch(type) {
    case MBEDGE:
    case MBTRI: {
      FENumeredDofMoFEMEntity_multiIndex::index<Composite_mi_tag>::type::iterator eiit;
      eiit = col_multiIndex->get<Composite_mi_tag>().find(boost::make_tuple(field_name,type,side_number));
      if(eiit == col_multiIndex->get<Composite_mi_tag>().end()) PetscFunctionReturn(0);//SETERRQ(PETSC_COMM_SELF,1,"no such ent");
      switch(type) {
	case MBEDGE: {
	  GlobIndices_EntType::iterator miit = col_edgesGlobIndices.find(eiit->get_MoFEMEntity_ptr());
	  if(miit == col_edgesGlobIndices.end()) SETERRQ(PETSC_COMM_SELF,1,"no such ent in FE");
	  ColGlobDofs = miit->second;
	} break;
	case MBTRI: {
	  GlobIndices_EntType::iterator miit = col_facesGlobIndices.find(eiit->get_MoFEMEntity_ptr());
	  if(miit == col_facesGlobIndices.end()) SETERRQ(PETSC_COMM_SELF,1,"no such ent in FE");
	  ColGlobDofs = miit->second;
	} break;
	default:
	SETERRQ(PETSC_COMM_SELF,1,"no implemented");
      }
    } break;
    case MBTET:
    case MBPRISM: {
      FENumeredDofMoFEMEntity_multiIndex::index<MoABEnt_mi_tag>::type::iterator eiit;
      eiit = col_multiIndex->get<MoABEnt_mi_tag>().find(fe_ent_ptr->get_ent());
      if(eiit == col_multiIndex->get<MoABEnt_mi_tag>().end()) PetscFunctionReturn(0);//SETERRQ(PETSC_COMM_SELF,1,"no such ent");
      GlobIndices_EntType::iterator miit = col_elemGlobIndices.find(eiit->get_MoFEMEntity_ptr());
      if(miit == col_elemGlobIndices.end()) SETERRQ(PETSC_COMM_SELF,1,"no such ent in FE");
      ColGlobDofs = miit->second;
    } break;
    default:
      SETERRQ(PETSC_COMM_SELF,1,"no implemented");
  }
  PetscFunctionReturn(0);
}
PetscErrorCode FEMethod_UpLevelStudent::GetDataVector(const string &field_name,ublas::vector<FieldData> &Data) {
  PetscFunctionBegin;
  MoFEMField_multiIndex::index<FieldName_mi_tag>::type::iterator fiit = moabfields->get<FieldName_mi_tag>().find(field_name);
  if(fiit==moabfields->get<FieldName_mi_tag>().end()) SETERRQ(PETSC_COMM_SELF,1,"no such field");
  Data_Type::iterator miit = data_nodes.find(fiit->get_MoFEMField_ptr());
  if(miit == data_nodes.end()) SETERRQ(PETSC_COMM_SELF,1,"no such field in FE");
  Data = miit->second;
  PetscFunctionReturn(0);
}
PetscErrorCode FEMethod_UpLevelStudent::GetDataVector(const string &field_name,EntityType type,ublas::vector<FieldData> &Data,int side_number) {
  PetscFunctionBegin;
  switch(type) {
    case MBEDGE:
    case MBTRI: {
      FENumeredDofMoFEMEntity_multiIndex::index<Composite_mi_tag>::type::iterator eiit;
      eiit = col_multiIndex->get<Composite_mi_tag>().find(boost::make_tuple(field_name,type,side_number));
      if(eiit == col_multiIndex->get<Composite_mi_tag>().end()) SETERRQ(PETSC_COMM_SELF,1,"no such ent");
      switch(type) {
	case MBEDGE: {
	  Data_EntType::iterator miit = data_edges.find(eiit->get_MoFEMEntity_ptr());
	  if(miit == data_edges.end()) SETERRQ(PETSC_COMM_SELF,1,"no such ent in FE");
	  Data = miit->second;
	} break;
	case MBTRI: {
	  Data_EntType::iterator miit = data_faces.find(eiit->get_MoFEMEntity_ptr());
	  if(miit == data_faces.end()) SETERRQ(PETSC_COMM_SELF,1,"no such ent in FE");
	  Data = miit->second;
	} break;
	default:
	SETERRQ(PETSC_COMM_SELF,1,"no implemented");
      }
    } break;
    case MBTET:
    case MBPRISM: {
      FENumeredDofMoFEMEntity_multiIndex::index<MoABEnt_mi_tag>::type::iterator eiit;
      eiit = col_multiIndex->get<MoABEnt_mi_tag>().find(fe_ent_ptr->get_ent());
      if(eiit == col_multiIndex->get<MoABEnt_mi_tag>().end()) SETERRQ(PETSC_COMM_SELF,1,"no such ent");
      Data_EntType::iterator miit = data_elem.find(eiit->get_MoFEMEntity_ptr());
      if(miit == data_elem.end()) SETERRQ(PETSC_COMM_SELF,1,"no such ent in FE");
      Data = miit->second;
    } break;
    default:
      SETERRQ(PETSC_COMM_SELF,1,"no implemented");
  }
  PetscFunctionReturn(0);
}
PetscErrorCode FEMethod_UpLevelStudent::GetGaussDataVector(const string &field_name,vector< ublas::vector<FieldData> > &Data) {
  PetscFunctionBegin;
  Data_at_Gauss_pt::iterator miit = data_at_gauss_pt.find(field_name);
  if(miit == data_at_gauss_pt.end()) SETERRQ(PETSC_COMM_SELF,1,"no such field in FE");
  Data = miit->second;
  PetscFunctionReturn(0);
}
PetscErrorCode FEMethod_UpLevelStudent::GetGaussDiffDataVector(const string &field_name,vector< ublas::matrix<FieldData> > &Data) {
  PetscFunctionBegin;
  DiffData_at_Gauss_pt::iterator miit = diff_data_at_gauss_pt.find(field_name);
  if(miit == diff_data_at_gauss_pt.end()) SETERRQ(PETSC_COMM_SELF,1,"no such field in FE");
  Data = miit->second;
  PetscFunctionReturn(0);
}
PetscErrorCode FEMethod_UpLevelStudent::GetGaussRowNMatrix(const string &field_name,vector< ublas::matrix<FieldData> > &NMatrix) {
  PetscFunctionBegin;
  MoFEMField_multiIndex::index<FieldName_mi_tag>::type::iterator fiit = moabfields->get<FieldName_mi_tag>().find(field_name);
  if(fiit==moabfields->get<FieldName_mi_tag>().end()) SETERRQ(PETSC_COMM_SELF,1,"no such field");
  N_Matrix_Type::iterator miit = row_N_Matrix_nodes.find(fiit->get_MoFEMField_ptr());
  if(miit == row_N_Matrix_nodes.end()) SETERRQ(PETSC_COMM_SELF,1,"no such field in FE");
  NMatrix = miit->second;
  PetscFunctionReturn(0);
}
PetscErrorCode FEMethod_UpLevelStudent::GetGaussRowNMatrix(const string &field_name,EntityType type,vector< ublas::matrix<FieldData> > &NMatrix,int side_number) {
  switch(type) {
    case MBEDGE:
    case MBTRI: {
      FENumeredDofMoFEMEntity_multiIndex::index<Composite_mi_tag>::type::iterator eiit;
      eiit = row_multiIndex->get<Composite_mi_tag>().find(boost::make_tuple(field_name,type,side_number));
      if(eiit == row_multiIndex->get<Composite_mi_tag>().end()) PetscFunctionReturn(0); //SETERRQ(PETSC_COMM_SELF,1,"no such ent");
      switch(type) {
	case MBEDGE: {
	  N_Matrix_EntType::iterator miit = row_N_Matrix_edges.find(eiit->get_MoFEMEntity_ptr());
	  if(miit == row_N_Matrix_edges.end()) SETERRQ(PETSC_COMM_SELF,1,"no such ent in FE");
	  NMatrix = miit->second;
	} break;
	case MBTRI: {
	  N_Matrix_EntType::iterator miit = row_N_Matrix_faces.find(eiit->get_MoFEMEntity_ptr());
	  if(miit == row_N_Matrix_faces.end()) SETERRQ(PETSC_COMM_SELF,1,"no such ent in FE");
	  NMatrix = miit->second;
	} break;
	default:
	SETERRQ(PETSC_COMM_SELF,1,"no implemented");
      }
    } break;
    case MBTET:
    case MBPRISM: {
      FENumeredDofMoFEMEntity_multiIndex::index<MoABEnt_mi_tag>::type::iterator eiit;
      eiit = row_multiIndex->get<MoABEnt_mi_tag>().find(fe_ent_ptr->get_ent());
      if(eiit == row_multiIndex->get<MoABEnt_mi_tag>().end()) PetscFunctionReturn(0);//SETERRQ(PETSC_COMM_SELF,1,"no such ent");
      N_Matrix_EntType::iterator miit = row_N_Matrix_elem.find(eiit->get_MoFEMEntity_ptr());
      if(miit == row_N_Matrix_elem.end()) SETERRQ(PETSC_COMM_SELF,1,"no such ent in FE");
      NMatrix = miit->second;
    } break;
    default:
      SETERRQ(PETSC_COMM_SELF,1,"no implemented");
  }
  PetscFunctionReturn(0);
}
PetscErrorCode FEMethod_UpLevelStudent::GetGaussColNMatrix(const string &field_name,vector< ublas::matrix<FieldData> > &NMatrix) {
  PetscFunctionBegin;
  MoFEMField_multiIndex::index<FieldName_mi_tag>::type::iterator fiit = moabfields->get<FieldName_mi_tag>().find(field_name);
  if(fiit==moabfields->get<FieldName_mi_tag>().end()) SETERRQ(PETSC_COMM_SELF,1,"no such field");
  N_Matrix_Type::iterator miit = col_N_Matrix_nodes.find(fiit->get_MoFEMField_ptr());
  if(miit == col_N_Matrix_nodes.end()) SETERRQ(PETSC_COMM_SELF,1,"no such field in FE");
  NMatrix = miit->second;
  PetscFunctionReturn(0);
}
PetscErrorCode FEMethod_UpLevelStudent::GetGaussColNMatrix(const string &field_name,EntityType type,vector< ublas::matrix<FieldData> > &NMatrix,int side_number) {
  switch(type) {
    case MBEDGE:
    case MBTRI: {
      FENumeredDofMoFEMEntity_multiIndex::index<Composite_mi_tag>::type::iterator eiit;
      eiit = col_multiIndex->get<Composite_mi_tag>().find(boost::make_tuple(field_name,type,side_number));
      if(eiit == col_multiIndex->get<Composite_mi_tag>().end()) PetscFunctionReturn(0);//SETERRQ(PETSC_COMM_SELF,1,"no such ent");
      switch(type) {
	case MBEDGE: {
	  N_Matrix_EntType::iterator miit = col_N_Matrix_edges.find(eiit->get_MoFEMEntity_ptr());
	  if(miit == col_N_Matrix_edges.end()) SETERRQ(PETSC_COMM_SELF,1,"no such ent in FE");
	  NMatrix = miit->second;
	} break;
	case MBTRI: {
	  N_Matrix_EntType::iterator miit = col_N_Matrix_faces.find(eiit->get_MoFEMEntity_ptr());
	  if(miit == col_N_Matrix_faces.end()) SETERRQ(PETSC_COMM_SELF,1,"no such ent in FE");
	  NMatrix = miit->second;
	} break;
	default:
	SETERRQ(PETSC_COMM_SELF,1,"no implemented");
      }
    } break;
    case MBTET:
    case MBPRISM: {
      FENumeredDofMoFEMEntity_multiIndex::index<MoABEnt_mi_tag>::type::iterator eiit;
      eiit = col_multiIndex->get<MoABEnt_mi_tag>().find(fe_ent_ptr->get_ent());
      if(eiit == col_multiIndex->get<MoABEnt_mi_tag>().end()) PetscFunctionReturn(0);//SETERRQ(PETSC_COMM_SELF,1,"no such ent");
      N_Matrix_EntType::iterator miit = col_N_Matrix_elem.find(eiit->get_MoFEMEntity_ptr());
      if(miit == col_N_Matrix_elem.end()) SETERRQ(PETSC_COMM_SELF,1,"no such ent in FE");
      NMatrix = miit->second;
    } break;
    default:
      SETERRQ(PETSC_COMM_SELF,1,"no implemented");
  }
  PetscFunctionReturn(0);
}
PetscErrorCode FEMethod_UpLevelStudent::GetGaussRowDiffNMatrix(const string &field_name,vector< ublas::matrix<FieldData> > &diffNMatrix) {
  PetscFunctionBegin;
  MoFEMField_multiIndex::index<FieldName_mi_tag>::type::iterator fiit = moabfields->get<FieldName_mi_tag>().find(field_name);
  if(fiit==moabfields->get<FieldName_mi_tag>().end()) SETERRQ(PETSC_COMM_SELF,1,"no such field");
  N_Matrix_Type::iterator miit = row_diffN_Matrix_nodes.find(fiit->get_MoFEMField_ptr());
  if(miit == row_N_Matrix_nodes.end()) SETERRQ(PETSC_COMM_SELF,1,"no such field in FE");
  diffNMatrix = miit->second;
  PetscFunctionReturn(0);
}
PetscErrorCode FEMethod_UpLevelStudent::GetGaussRowDiffNMatrix(const string &field_name,EntityType type,vector< ublas::matrix<FieldData> > &diffNMatrix,int side_number) {
  switch(type) {
    case MBEDGE:
    case MBTRI: {
      FENumeredDofMoFEMEntity_multiIndex::index<Composite_mi_tag>::type::iterator eiit;
      eiit = row_multiIndex->get<Composite_mi_tag>().find(boost::make_tuple(field_name,type,side_number));
      if(eiit == row_multiIndex->get<Composite_mi_tag>().end()) SETERRQ(PETSC_COMM_SELF,1,"no such ent");
      switch(type) {
	case MBEDGE: {
	  N_Matrix_EntType::iterator miit = row_diffN_Matrix_edges.find(eiit->get_MoFEMEntity_ptr());
	  if(miit == row_N_Matrix_edges.end()) SETERRQ(PETSC_COMM_SELF,1,"no such ent in FE");
	  diffNMatrix = miit->second;
	} break;
	case MBTRI: {
	  N_Matrix_EntType::iterator miit = row_diffN_Matrix_faces.find(eiit->get_MoFEMEntity_ptr());
	  if(miit == row_N_Matrix_faces.end()) SETERRQ(PETSC_COMM_SELF,1,"no such ent in FE");
	  diffNMatrix = miit->second;
	} break;
	default:
	SETERRQ(PETSC_COMM_SELF,1,"no implemented");
      }
    } break;
    case MBTET:
    case MBPRISM: {
      FENumeredDofMoFEMEntity_multiIndex::index<MoABEnt_mi_tag>::type::iterator eiit;
      eiit = row_multiIndex->get<MoABEnt_mi_tag>().find(fe_ent_ptr->get_ent());
      if(eiit == row_multiIndex->get<MoABEnt_mi_tag>().end()) SETERRQ(PETSC_COMM_SELF,1,"no such ent");
      N_Matrix_EntType::iterator miit = row_diffN_Matrix_elem.find(eiit->get_MoFEMEntity_ptr());
      if(miit == row_diffN_Matrix_elem.end()) SETERRQ(PETSC_COMM_SELF,1,"no such ent in FE");
      diffNMatrix = miit->second;
    } break;
    default:
      SETERRQ(PETSC_COMM_SELF,1,"no implemented");
  }
  PetscFunctionReturn(0);
}
PetscErrorCode FEMethod_UpLevelStudent::GetGaussColDiffNMatrix(const string &field_name,vector< ublas::matrix<FieldData> > &diffNMatrix) {
  PetscFunctionBegin;
  MoFEMField_multiIndex::index<FieldName_mi_tag>::type::iterator fiit = moabfields->get<FieldName_mi_tag>().find(field_name);
  if(fiit==moabfields->get<FieldName_mi_tag>().end()) SETERRQ(PETSC_COMM_SELF,1,"no such field");
  N_Matrix_Type::iterator miit = col_diffN_Matrix_nodes.find(fiit->get_MoFEMField_ptr());
  if(miit == col_N_Matrix_nodes.end()) SETERRQ(PETSC_COMM_SELF,1,"no such field in FE");
  diffNMatrix = miit->second;
  PetscFunctionReturn(0);
}
PetscErrorCode FEMethod_UpLevelStudent::GetGaussColDiffNMatrix(const string &field_name,EntityType type,vector< ublas::matrix<FieldData> > &diffNMatrix,int side_number) {
  switch(type) {
    case MBEDGE:
    case MBTRI: {
      FENumeredDofMoFEMEntity_multiIndex::index<Composite_mi_tag>::type::iterator eiit;
      eiit = col_multiIndex->get<Composite_mi_tag>().find(boost::make_tuple(field_name,type,side_number));
      if(eiit == col_multiIndex->get<Composite_mi_tag>().end()) SETERRQ(PETSC_COMM_SELF,1,"no such ent");
      switch(type) {
	case MBEDGE: {
	  N_Matrix_EntType::iterator miit = col_diffN_Matrix_edges.find(eiit->get_MoFEMEntity_ptr());
	  if(miit == col_N_Matrix_edges.end()) SETERRQ(PETSC_COMM_SELF,1,"no such ent in FE");
	  diffNMatrix = miit->second;
	} break;
	case MBTRI: {
	  N_Matrix_EntType::iterator miit = col_diffN_Matrix_faces.find(eiit->get_MoFEMEntity_ptr());
	  if(miit == col_N_Matrix_faces.end()) SETERRQ(PETSC_COMM_SELF,1,"no such ent in FE");
	  diffNMatrix = miit->second;
	} break;
	default:
	SETERRQ(PETSC_COMM_SELF,1,"no implemented");
      }
    } break;
    case MBTET:
    case MBPRISM: {
      FENumeredDofMoFEMEntity_multiIndex::index<MoABEnt_mi_tag>::type::iterator eiit;
      eiit = col_multiIndex->get<MoABEnt_mi_tag>().find(fe_ent_ptr->get_ent());
      if(eiit == col_multiIndex->get<MoABEnt_mi_tag>().end()) SETERRQ(PETSC_COMM_SELF,1,"no such ent");
      N_Matrix_EntType::iterator miit = col_diffN_Matrix_elem.find(eiit->get_MoFEMEntity_ptr());
      if(miit == col_N_Matrix_elem.end()) SETERRQ(PETSC_COMM_SELF,1,"no such ent in FE");
      diffNMatrix = miit->second;
    } break;
    default:
      SETERRQ(PETSC_COMM_SELF,1,"no implemented");
  }
  PetscFunctionReturn(0);
}
PetscErrorCode FEMethod_UpLevelStudent::MakeBMatrix3D(
  const string &field_name,vector<ublas::matrix<FieldData> > &diffNMatrix,vector<ublas::matrix<FieldData> > &BMatrix) {
  PetscFunctionBegin;
  MoFEMField_multiIndex::index<FieldName_mi_tag>::type::iterator fiit = moabfields->get<FieldName_mi_tag>().find(field_name);
  if(fiit==moabfields->get<FieldName_mi_tag>().end()) SETERRQ(PETSC_COMM_SELF,1,"no such field");
  if(fiit->get_space() != H1) SETERRQ(PETSC_COMM_SELF,1,"it has to be H1 space");
  if(fiit->get_max_rank() != 3) SETERRQ(PETSC_COMM_SELF,1,"it has to be rank 3");
  const int g_dim = get_dim_gNTET();
  BMatrix.resize(g_dim);
  for(int gg = 0;gg<g_dim;gg++) {
    ublas::matrix<FieldData> &diffMat = diffNMatrix[gg];
    ublas::matrix<FieldData> &BMat = BMatrix[gg];
    int m = diffMat.size1();
    int n = diffMat.size2();
    if(m != 9)  SETERRQ(PETSC_COMM_SELF,1,"wrong matrix size");
    BMat.resize(6,n);
    // page 30 CHAPTER 6. DISPLACEMENT METHODS, FEAP Version 7.3 Theory Manual Robert L. Taylor
    ublas::matrix_row<ublas::matrix<FieldData> >(BMat,0) = ublas::matrix_row<ublas::matrix<FieldData> >(diffMat,0); //dX/dx
    ublas::matrix_row<ublas::matrix<FieldData> >(BMat,1) = ublas::matrix_row<ublas::matrix<FieldData> >(diffMat,4); //dY/dy
    ublas::matrix_row<ublas::matrix<FieldData> >(BMat,2) = ublas::matrix_row<ublas::matrix<FieldData> >(diffMat,8); //dZ/dz
    ublas::matrix_row<ublas::matrix<FieldData> >(BMat,3) = //dX/dy+dY/dx
      ublas::matrix_row<ublas::matrix<FieldData> >(diffMat,1)+ublas::matrix_row<ublas::matrix<FieldData> >(diffMat,3);
    ublas::matrix_row<ublas::matrix<FieldData> >(BMat,4) = //dY/dz+dZ/dy
      ublas::matrix_row<ublas::matrix<FieldData> >(diffMat,5)+ublas::matrix_row<ublas::matrix<FieldData> >(diffMat,7);
    ublas::matrix_row<ublas::matrix<FieldData> >(BMat,5) = //dX/dz+dZ/dx
      ublas::matrix_row<ublas::matrix<FieldData> >(diffMat,2)+ublas::matrix_row<ublas::matrix<FieldData> >(diffMat,6);
  }
  PetscFunctionReturn(0);
}
PetscErrorCode FEMethod_UpLevelStudent::GetGaussRowFaceNMatrix(
  EntityHandle ent,const string &field_name,
  vector< ublas::matrix<FieldData> > &NMatrix,
  EntityType type,EntityHandle edge_handle) {
  PetscFunctionBegin;
  N_Matrix_Type N_Matrix_nodes;
  N_Matrix_EntType N_Matrix_edges;
  N_Matrix_EntType N_Matrix_faces;
  ierr = GetNMatrix_at_FaceGaussPoint(ent,field_name,
    row_nodesGlobIndices,row_edgesGlobIndices,row_facesGlobIndices,
    N_Matrix_nodes,N_Matrix_edges,N_Matrix_faces,
    type,edge_handle); CHKERRQ(ierr);
  MoFEMField_multiIndex::index<FieldName_mi_tag>::type::iterator fiit = moabfields->get<FieldName_mi_tag>().find(field_name);
  if(fiit==moabfields->get<FieldName_mi_tag>().end()) SETERRQ(PETSC_COMM_SELF,1,"no such field");
  switch (type) {
    case MBVERTEX: {
      N_Matrix_Type::iterator miit = N_Matrix_nodes.find(fiit->get_MoFEMField_ptr());
      if(miit == col_N_Matrix_nodes.end()) SETERRQ(PETSC_COMM_SELF,1,"no such field in FE");
      NMatrix = miit->second;
      }
      break;
    case MBTRI: {
      int side_number;
      try {
	side_number = fe_ent_ptr->get_side_number_ptr(moab,ent)->side_number;
      } catch (const char* msg) {
	  SETERRQ(PETSC_COMM_SELF,1,msg);
      }
      FENumeredDofMoFEMEntity_multiIndex::index<Composite_mi_tag>::type::iterator fiiit;
      fiiit = row_multiIndex->get<Composite_mi_tag>().find(boost::make_tuple(field_name,MBTRI,side_number));
      if(fiiit == row_multiIndex->get<Composite_mi_tag>().end()) SETERRQ1(PETSC_COMM_SELF,1,"no such ent (side_number = %u)",side_number);
      N_Matrix_EntType::iterator miit = N_Matrix_faces.find(fiiit->get_MoFEMEntity_ptr());
      if(miit == N_Matrix_faces.end()) SETERRQ(PETSC_COMM_SELF,1,"no such field in FE");
      NMatrix = miit->second;
      }
      break;
    case MBEDGE: {
      int side_number;
      try {
	side_number = fe_ent_ptr->get_side_number_ptr(moab,edge_handle)->side_number;
      } catch (const char* msg) {
	  SETERRQ(PETSC_COMM_SELF,1,msg);
      }
      FENumeredDofMoFEMEntity_multiIndex::index<Composite_mi_tag>::type::iterator eiiit;
      eiiit = row_multiIndex->get<Composite_mi_tag>().find(boost::make_tuple(field_name,MBEDGE,side_number));
      if(eiiit == row_multiIndex->get<Composite_mi_tag>().end()) SETERRQ1(PETSC_COMM_SELF,1,"no such ent (side_number = %u)",side_number);
      N_Matrix_EntType::iterator miit = N_Matrix_edges.find(eiiit->get_MoFEMEntity_ptr());
      if(miit == N_Matrix_faces.end()) SETERRQ(PETSC_COMM_SELF,1,"no such field in FE");
      NMatrix = miit->second;
      }
      break;
    default:
      SETERRQ(PETSC_COMM_SELF,1,"no implemented");
  }
  
  PetscFunctionReturn(0);
}


}

