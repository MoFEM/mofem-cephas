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

#ifndef __MOABSURFACECONSTRAINS_HPP__
#define __MOABSURFACECONSTRAINS_HPP__

namespace MoFEM {

struct C_SURFACE_FEMethod:public moabField::FEMethod {
  ErrorCode rval;
  PetscErrorCode ierr;
  Interface& moab;

  Mat C;
  Range skin_faces;
  Tag th_material_normal;
  vector<double> diffNTRI;
  vector<double> g_NTRI3;
  const double *G_TRI_W;
  
  void run_in_constructor() {
    diffNTRI.resize(6);
    ShapeDiffMBTRI(&diffNTRI[0]);
    g_NTRI3.resize(3*3);
    ShapeMBTRI(&g_NTRI3[0],G_TRI_X3,G_TRI_Y3,3);
    G_TRI_W = G_TRI_W3;
    double def_VAL[3*9];
    fill(&def_VAL[0],&def_VAL[3*9],0);
    rval = moab.tag_get_handle("MATERIAL_NORMAL",3,MB_TYPE_DOUBLE,th_material_normal,MB_TAG_CREAT|MB_TAG_SPARSE,&def_VAL); CHKERR_THROW(rval);
  }
  
  C_SURFACE_FEMethod(Interface& _moab,EntityHandle skin_faces_meshset,Mat _C,int _verbose = 0): 
    FEMethod(),moab(_moab),C(_C) {
    run_in_constructor();
    rval = moab.get_entities_by_type(skin_faces_meshset,MBTRI,skin_faces,true);  CHKERR_THROW(rval);
  }
  C_SURFACE_FEMethod(Interface& _moab,Range &_skin_faces,Mat _C,int _verbose = 0): 
    FEMethod(),moab(_moab),C(_C),skin_faces(_skin_faces) {
    run_in_constructor();
  }
  C_SURFACE_FEMethod(Interface& _moab,Mat _C,int _verbose = 0): 
    FEMethod(),moab(_moab),C(_C) {} 

  PetscErrorCode preProcess() {
    PetscFunctionBegin;
    PetscFunctionReturn(0);
  }

  ublas::matrix<double> C_MAT_ELEM;
  ublas::vector<DofIdx> ent_global_col_indices,ent_global_row_indices;
  ublas::vector<double,ublas::bounded_array<double,9> > ent_dofs_data;
  ublas::vector<double,ublas::bounded_array<double,3> > ent_normal_map;
  ublas::vector<double,ublas::bounded_array<double,3> > ent_normal_map0;
  ublas::vector<double> coords;
  virtual PetscErrorCode Integrate() {
    PetscFunctionBegin;
    C_MAT_ELEM.resize(3,9);
    ublas::noalias(C_MAT_ELEM) = ublas::zero_matrix<double>(3,9);
    double area0 = norm_2(ent_normal_map0);
    double area = norm_2(ent_normal_map);
    for(int gg = 0;gg<g_NTRI3.size()/3;gg++) {
	for(int nn = 0;nn<3;nn++) {
	  for(int dd = 0;dd<3;dd++) {
	    for(int nnn = 0;nnn<3;nnn++) {
	      C_MAT_ELEM(nn,3*nnn+dd) += G_TRI_W[gg]*ent_normal_map[dd]*g_NTRI3[3*gg+nn]*g_NTRI3[3*gg+nnn]*(area0/area);
	    }
	  }
	}
    }
    /*cerr << "ROWS " << ent_global_row_indices << endl;
    cerr << "COLS " << ent_global_col_indices << endl;
    cerr << "NORMAL " << ent_normal_map << endl;
    cerr << "MAT " << C_MAT_ELEM << endl;*/
    ierr = MatSetValues(C,
      ent_global_row_indices.size(),&(ent_global_row_indices.data()[0]),
      ent_global_col_indices.size(),&(ent_global_col_indices.data()[0]),
      &(C_MAT_ELEM.data())[0],ADD_VALUES); CHKERRQ(ierr);
    PetscFunctionReturn(0);
  }

  PetscErrorCode operator()() {
    PetscFunctionBegin;
    SideNumber_multiIndex &side_table = fe_ptr->get_side_number_table();
    SideNumber_multiIndex::nth_index<1>::type::iterator siit = side_table.get<1>().lower_bound(boost::make_tuple(MBTRI,0));
    SideNumber_multiIndex::nth_index<1>::type::iterator hi_siit = side_table.get<1>().upper_bound(boost::make_tuple(MBTRI,4));
    for(;siit!=hi_siit;siit++) {
      EntityHandle face = siit->ent;
      if(find(skin_faces.begin(),skin_faces.end(),face)==skin_faces.end()) continue;
      ent_dofs_data.resize(9);
      ent_global_col_indices.resize(9);
      ent_global_row_indices.resize(3);
      const EntityHandle* conn_face; 
      int num_nodes; 
      rval = moab.get_connectivity(face,conn_face,num_nodes,true); CHKERR_PETSC(rval);
      for(int nn = 0;nn<num_nodes;nn++) {
	FENumeredDofMoFEMEntity_multiIndex::index<Composite_mi_tag3>::type::iterator dit,hi_dit;
	dit = row_multiIndex->get<Composite_mi_tag3>().lower_bound(boost::make_tuple("LAMBDA_SURFACE",conn_face[nn]));
	hi_dit = row_multiIndex->get<Composite_mi_tag3>().upper_bound(boost::make_tuple("LAMBDA_SURFACE",conn_face[nn]));
	if(distance(dit,hi_dit)==0) {
	  ent_global_row_indices[nn] = -1;
	} else {
	  if(distance(dit,hi_dit)!=1) SETERRQ(PETSC_COMM_SELF,1,"data inconsistency, number of dof on node for LAMBDA_SURFACE should be 1");
	  int global_idx = dit->get_petsc_gloabl_dof_idx();
	  ent_global_row_indices[nn] = global_idx;
	  int local_idx = dit->get_petsc_local_dof_idx();
	  if(local_idx<0) SETERRQ(PETSC_COMM_SELF,1,"data inconsistency, negative index of local dofs on element");
	}
	dit = col_multiIndex->get<Composite_mi_tag3>().lower_bound(boost::make_tuple("MESH_NODE_POSITIONS",conn_face[nn]));
	hi_dit = col_multiIndex->get<Composite_mi_tag3>().upper_bound(boost::make_tuple("MESH_NODE_POSITIONS",conn_face[nn]));
	if(distance(dit,hi_dit)!=3) SETERRQ(PETSC_COMM_SELF,1,"data inconsistency, number of dof on node for MESH_NODE_POSITIONS should be 3");
	  for(;dit!=hi_dit;dit++) {
	    int global_idx = dit->get_petsc_gloabl_dof_idx();
	    assert(nn*3+dit->get_dof_rank()<9);
	    ent_global_col_indices[nn*3+dit->get_dof_rank()] = global_idx;
	    int local_idx = dit->get_petsc_local_dof_idx();
	    if(local_idx<0) SETERRQ(PETSC_COMM_SELF,1,"data inconsistency, negative index of local dofs on element");
	    ent_dofs_data[nn*3+dit->get_dof_rank()] = dit->get_FieldData();
	}
      }
      coords.resize(num_nodes*3);
      rval = moab.get_coords(conn_face,num_nodes,&*coords.data().begin()); CHKERR_PETSC(rval);
      ent_normal_map0.resize(3);
      ierr = ShapeFaceNormalMBTRI(&diffNTRI[0],&coords.data()[0],&ent_normal_map0.data()[0]); CHKERRQ(ierr);
      ent_normal_map.resize(3);
      ierr = ShapeFaceNormalMBTRI(&diffNTRI[0],&ent_dofs_data.data()[0],&ent_normal_map.data()[0]); CHKERRQ(ierr);
      ent_normal_map *= siit->sense;
      rval = moab.tag_set_data(th_material_normal,&face,1,&ent_normal_map.data()[0]); CHKERR_PETSC(rval);
      ierr = this->Integrate(); CHKERRQ(ierr);
    }
    PetscFunctionReturn(0);
  }
  PetscErrorCode postProcess() {
    PetscFunctionBegin;
    PetscFunctionReturn(0);
  }

};

struct g_SURFACE_FEMethod: public C_SURFACE_FEMethod {
  
  Vec g;
  g_SURFACE_FEMethod(Interface& _moab,EntityHandle skin_faces_meshset,Vec _g,int _verbose = 0): 
    C_SURFACE_FEMethod(_moab,skin_faces,PETSC_NULL,_verbose),g(_g) {
    VecSetOption(g, VEC_IGNORE_NEGATIVE_INDICES,PETSC_TRUE); 
  }
  g_SURFACE_FEMethod(Interface& _moab,Range &_skin_faces,Vec _g,int _verbose = 0): 
    C_SURFACE_FEMethod(_moab,_skin_faces,PETSC_NULL,_verbose),g(_g) {
    VecSetOption(g, VEC_IGNORE_NEGATIVE_INDICES,PETSC_TRUE); 
  };

  ublas::vector<double,ublas::bounded_array<double,3> > g_VEC_ELEM;
  PetscErrorCode Integrate() {
    PetscFunctionBegin;
    g_VEC_ELEM.resize(3);
    ublas::noalias(g_VEC_ELEM) = ublas::zero_vector<double>(3);
    double area0 = norm_2(ent_normal_map0);
    double area = norm_2(ent_normal_map);
    for(int gg = 0;gg<g_NTRI3.size()/3;gg++) {
	for(int nn = 0;nn<3;nn++) {
	  for(int dd = 0;dd<3;dd++) {
	    double X0_dd = cblas_ddot(3,&g_NTRI3[0],1,&coords.data()[dd],3);
	    double X_dd = cblas_ddot(3,&g_NTRI3[0],1,&ent_dofs_data.data()[dd],3);
	    g_VEC_ELEM[nn] += G_TRI_W[gg]*g_NTRI3[3*gg+nn]*ent_normal_map[dd]*(X0_dd - X_dd)*(area0/area);
	  }
	}
    }
    if(ent_global_row_indices.size()!=g_VEC_ELEM.size()) SETERRQ(PETSC_COMM_SELF,1,"data inconsistency");
    ierr = VecSetValues(g,
      ent_global_row_indices.size(),&(ent_global_row_indices.data()[0]),
      &(g_VEC_ELEM.data())[0],ADD_VALUES); CHKERRQ(ierr);
    PetscFunctionReturn(0);
  }

};

struct C_EDGE_FEMethod:public C_SURFACE_FEMethod {

  Range edges,faces;
  C_EDGE_FEMethod(Interface& _moab,EntityHandle edges_meshset,Mat _C,int _verbose = 0): 
    C_SURFACE_FEMethod(_moab,_C,_verbose) {
    run_in_constructor();
    rval = moab.get_entities_by_type(edges_meshset,MBEDGE,edges,true);  CHKERR_THROW(rval);
    rval = moab.get_entities_by_type(edges_meshset,MBTRI,faces,true);  CHKERR_THROW(rval);
  }
  C_EDGE_FEMethod(Interface& _moab,Range &_faces,Range &_edges,Mat _C,int _verbose = 0): 
  C_SURFACE_FEMethod(_moab,_C,_verbose),faces(_faces),edges(_edges) {
    run_in_constructor();
  }

  map<EntityHandle,vector<EntityHandle> > edges_face_map;
  PetscErrorCode preProcess() {
    PetscFunctionBegin;
    edges_face_map.clear();
    PetscFunctionReturn(0);
  }

  PetscErrorCode operator()() {
    PetscFunctionBegin;
    SideNumber_multiIndex &side_table = fe_ptr->get_side_number_table();
    SideNumber_multiIndex::nth_index<1>::type::iterator siit = side_table.get<1>().lower_bound(boost::make_tuple(MBTRI,0));
    SideNumber_multiIndex::nth_index<1>::type::iterator hi_siit = side_table.get<1>().upper_bound(boost::make_tuple(MBTRI,4));
    for(;siit!=hi_siit;siit++) {
      EntityHandle face = siit->ent;
      if(find(faces.begin(),faces.end(),face)==faces.end()) continue;
      Range adj_edges;
      rval = moab.get_adjacencies(&face,1,1,false,adj_edges,Interface::UNION); CHKERR_PETSC(rval);
      adj_edges = intersect(edges,adj_edges);
      Range::iterator eit = adj_edges.begin();
      for(;eit!=adj_edges.end();eit++) {
	map<EntityHandle,vector<EntityHandle> >::iterator mit = edges_face_map.find(*eit);
	int eq_nb = -1;
	if(mit == edges_face_map.end()) {
	  edges_face_map[*eit].push_back(face);
	  eq_nb = 0;
	} else {
	  if(mit->second.size()==1) {
	    if(mit->second[0]==face) {
	      eq_nb = 0;
	    } else {
	      mit->second.push_back(face);
	      eq_nb = 1;
	    }
	  } else {
	    if(mit->second[0]==face) eq_nb = 0;
	    else if(mit->second[1]==face) eq_nb = 1;
	    else SETERRQ(PETSC_COMM_SELF,1,"something is wrong");
	  }
	}
	if(eq_nb == -1) SETERRQ(PETSC_COMM_SELF,1,"something is wrong");
	ent_dofs_data.resize(9);
	ent_global_col_indices.resize(9);
	ent_global_row_indices.resize(3);
	const EntityHandle* conn_face; 
	int num_nodes; 
	rval = moab.get_connectivity(face,conn_face,num_nodes,true); CHKERR_PETSC(rval);
	for(int nn = 0;nn<num_nodes;nn++) {
	  FENumeredDofMoFEMEntity_multiIndex::index<Composite_mi_tag3>::type::iterator dit,hi_dit;
	  dit = row_multiIndex->get<Composite_mi_tag3>().lower_bound(boost::make_tuple("LAMBDA_EDGE",conn_face[nn]));
	  hi_dit = row_multiIndex->get<Composite_mi_tag3>().upper_bound(boost::make_tuple("LAMBDA_EDGE",conn_face[nn]));
	  if(distance(dit,hi_dit)==0) {
	    ent_global_row_indices[nn]=-1;
	  } else {
	    if(distance(dit,hi_dit)!=2) SETERRQ(PETSC_COMM_SELF,1,"data inconsistency, number of dof on node for LAMBDA_EDGE should be 2");
	    if(eq_nb==1) dit++;
	    int global_idx = dit->get_petsc_gloabl_dof_idx();
	    ent_global_row_indices[nn] = global_idx;
	    int local_idx = dit->get_petsc_local_dof_idx();
	    if(local_idx<0) SETERRQ(PETSC_COMM_SELF,1,"data inconsistency, negative index of local dofs on element");
	  }
	  dit = col_multiIndex->get<Composite_mi_tag3>().lower_bound(boost::make_tuple("MESH_NODE_POSITIONS",conn_face[nn]));
	  hi_dit = col_multiIndex->get<Composite_mi_tag3>().upper_bound(boost::make_tuple("MESH_NODE_POSITIONS",conn_face[nn]));
	  if(distance(dit,hi_dit)!=3) SETERRQ(PETSC_COMM_SELF,1,"data inconsistency, number of dof on node for MESH_NODE_POSITIONS should be 3");
	  for(;dit!=hi_dit;dit++) {
	      int global_idx = dit->get_petsc_gloabl_dof_idx();
	      assert(nn*3+dit->get_dof_rank()<9);
	      ent_global_col_indices[nn*3+dit->get_dof_rank()] = global_idx;
	      int local_idx = dit->get_petsc_local_dof_idx();
	      if(local_idx<0) SETERRQ(PETSC_COMM_SELF,1,"data inconsistency, negative index of local dofs on element");
	      ent_dofs_data[nn*3+dit->get_dof_rank()] = dit->get_FieldData();
	  }
	}
	coords.resize(num_nodes*3);
	rval = moab.get_coords(conn_face,num_nodes,&*coords.data().begin()); CHKERR_PETSC(rval);
	ent_normal_map0.resize(3);
	ierr = ShapeFaceNormalMBTRI(&diffNTRI[0],&coords.data()[0],&ent_normal_map0.data()[0]); CHKERRQ(ierr);
	ent_normal_map.resize(3);
	ierr = ShapeFaceNormalMBTRI(&diffNTRI[0],&ent_dofs_data.data()[0],&ent_normal_map.data()[0]); CHKERRQ(ierr);
	ent_normal_map *= siit->sense;
	ierr = this->Integrate(); CHKERRQ(ierr);
      }
    }
    PetscFunctionReturn(0);
  }

};

struct g_EDGE_FEMethod: public C_EDGE_FEMethod {
  
  Vec g;
  g_EDGE_FEMethod(Interface& _moab,EntityHandle edges_meshset,Vec _g,int _verbose = 0): 
    C_EDGE_FEMethod(_moab,edges_meshset,PETSC_NULL,_verbose),g(_g) {
    VecSetOption(g, VEC_IGNORE_NEGATIVE_INDICES,PETSC_TRUE); 
  }
  g_EDGE_FEMethod(Interface& _moab,Range &_faces,Range &_edges,Vec _g,int _verbose = 0): 
    C_EDGE_FEMethod(_moab,_faces,_edges,PETSC_NULL,_verbose),g(_g) {
    VecSetOption(g, VEC_IGNORE_NEGATIVE_INDICES,PETSC_TRUE); 
  };

  ublas::vector<double,ublas::bounded_array<double,3> > g_VEC_ELEM;
  PetscErrorCode Integrate() {
    PetscFunctionBegin;
    g_VEC_ELEM.resize(3);
    ublas::noalias(g_VEC_ELEM) = ublas::zero_vector<double>(3);
    double area0 = norm_2(ent_normal_map0);
    double area = norm_2(ent_normal_map);
    for(int gg = 0;gg<g_NTRI3.size()/3;gg++) {
      for(int nn = 0;nn<3;nn++) {
	for(int dd = 0;dd<3;dd++) {
	  double X0_dd = cblas_ddot(3,&g_NTRI3[0],1,&coords.data()[dd],3);
	  double X_dd = cblas_ddot(3,&g_NTRI3[0],1,&ent_dofs_data.data()[dd],3);
	  g_VEC_ELEM[nn] += G_TRI_W[gg]*g_NTRI3[3*gg+nn]*ent_normal_map[dd]*(X0_dd - X_dd)*(area0/area);
	}
      }
    }
    ierr = VecSetValues(g,
      ent_global_row_indices.size(),&(ent_global_row_indices.data()[0]),
      &(g_VEC_ELEM.data())[0],ADD_VALUES); CHKERRQ(ierr);
    PetscFunctionReturn(0);
  }

};

struct C_CORNER_FEMethod:public moabField::FEMethod {
  ErrorCode rval;
  PetscErrorCode ierr;
  Interface& moab;

  Mat C;
  Range corners;
  C_CORNER_FEMethod(Interface& _moab,Range _corners,Mat _C,int _verbose = 0): 
    FEMethod(),moab(_moab),C(_C),corners(_corners) {}

  PetscErrorCode preProcess() {
    PetscFunctionBegin;
    PetscFunctionReturn(0);
  }
  ublas::matrix<double> C_MAT_ELEM;
  vector<DofIdx> ent_global_col_indices,ent_global_row_indices;
  ublas::vector<double,ublas::bounded_array<double,9> > ent_dofs_data;
  ublas::vector<double,ublas::bounded_array<double,3> > ent_normal_map;
  ublas::vector<double,ublas::bounded_array<double,3> > coords; 
  virtual PetscErrorCode Integrate() {
    PetscFunctionBegin;
    C_MAT_ELEM.resize(3,3);
    for(int nn = 0;nn<3;nn++) {
	for(int dd = 0;dd<3;dd++) {
	  if(nn!=dd) C_MAT_ELEM(nn,dd) = 0;
	  else C_MAT_ELEM(nn,dd) = 1;
	}
    }
    //cerr << C_MAT_ELEM << endl;
    ierr = MatSetValues(C,
      ent_global_row_indices.size(),&(ent_global_row_indices[0]),
      ent_global_col_indices.size(),&(ent_global_col_indices[0]),
      &(C_MAT_ELEM.data())[0],INSERT_VALUES); CHKERRQ(ierr);
    PetscFunctionReturn(0);
  }
  PetscErrorCode operator()() {
    PetscFunctionBegin;
    SideNumber_multiIndex &side_table = fe_ptr->get_side_number_table();
    SideNumber_multiIndex::nth_index<1>::type::iterator siit = side_table.get<1>().lower_bound(boost::make_tuple(MBVERTEX,0));
    SideNumber_multiIndex::nth_index<1>::type::iterator hi_siit = side_table.get<1>().upper_bound(boost::make_tuple(MBVERTEX,4));
    for(;siit!=hi_siit;siit++) {
      EntityHandle node = siit->ent;
      if(corners.end()==find(corners.begin(),corners.end(),node)) continue;
      ent_dofs_data.resize(3);
      ent_global_col_indices.resize(3);
      ent_global_row_indices.resize(3);
      FENumeredDofMoFEMEntity_multiIndex::index<Composite_mi_tag3>::type::iterator dit,hi_dit;
      dit = col_multiIndex->get<Composite_mi_tag3>().lower_bound(boost::make_tuple("MESH_NODE_POSITIONS",node));
      hi_dit = col_multiIndex->get<Composite_mi_tag3>().upper_bound(boost::make_tuple("MESH_NODE_POSITIONS",node));
      if(distance(dit,hi_dit)!=3) SETERRQ(PETSC_COMM_SELF,1,"data inconsistency, number of dof on node for MESH_NODE_POSITIONS should be 3");
      for(;dit!=hi_dit;dit++) {
	int global_idx = dit->get_petsc_gloabl_dof_idx();
	ent_global_col_indices[dit->get_dof_rank()] = global_idx;
	int local_idx = dit->get_petsc_local_dof_idx();
	if(local_idx<0) SETERRQ(PETSC_COMM_SELF,1,"data inconsistency, negative index of local dofs on element");
	ent_dofs_data[dit->get_dof_rank()] = dit->get_FieldData();
      }
      dit = row_multiIndex->get<Composite_mi_tag3>().lower_bound(boost::make_tuple("LAMBDA_CORNER",node));
      hi_dit = row_multiIndex->get<Composite_mi_tag3>().upper_bound(boost::make_tuple("LAMBDA_CORNER",node));
      if(distance(dit,hi_dit)!=3) SETERRQ1(PETSC_COMM_SELF,1,"data inconsistency, number of dof on node for LAMBDA_CORNER should be 3, but is %d",distance(dit,hi_dit));
      for(;dit!=hi_dit;dit++) {
	int global_idx = dit->get_petsc_gloabl_dof_idx();
	ent_global_row_indices[dit->get_dof_rank()] = global_idx;
	int local_idx = dit->get_petsc_local_dof_idx();
	if(local_idx<0) SETERRQ(PETSC_COMM_SELF,1,"data inconsistency, negative index of local dofs on element");
      }
      coords.resize(3);
      rval = moab.get_coords(&node,1,&*coords.data().begin()); CHKERR_PETSC(rval);
      ierr = this->Integrate(); CHKERRQ(ierr);
    }
    PetscFunctionReturn(0);
  }
  PetscErrorCode postProcess() {
    PetscFunctionBegin;
    PetscFunctionReturn(0);
  }

};

struct g_CORNER_FEMethod: public C_CORNER_FEMethod {
  
  Vec g;
  g_CORNER_FEMethod(Interface& _moab,Range _corners,Vec _g,int _verbose = 0): 
    C_CORNER_FEMethod(_moab,_corners,PETSC_NULL,_verbose),g(_g) {
    VecSetOption(g, VEC_IGNORE_NEGATIVE_INDICES,PETSC_TRUE); 
  }

  ublas::vector<double,ublas::bounded_array<double,3> > g_VEC_ELEM;
  PetscErrorCode Integrate() {
    PetscFunctionBegin;
    g_VEC_ELEM.resize(3);
    for(int nn = 0;nn<3;nn++) {
      g_VEC_ELEM[nn] = coords[nn] - ent_dofs_data[nn];
    }
    ierr = VecSetValues(g,
      ent_global_row_indices.size(),&(ent_global_row_indices[0]),
      &(g_VEC_ELEM.data())[0],INSERT_VALUES); CHKERRQ(ierr);
    PetscFunctionReturn(0);
  }

};


}

#endif //__MOABSURFACECONSTRAINS_HPP__
