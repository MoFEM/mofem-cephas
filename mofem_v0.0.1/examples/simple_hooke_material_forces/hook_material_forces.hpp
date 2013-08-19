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

#ifndef __NONLINEAR_ELASTICITY_HPP__
#define __NONLINEAR_ELASTICITY_HPP__

#include "moabField.hpp"
#include "moabField_Core.hpp"
#include "moabFEMethod_UpLevelStudent.hpp"
#include "cholesky.hpp"
#include <petscksp.h>

#include "moabSnes.hpp"
#include "moabFEMethod_ComplexForLazy.hpp"
#include "moabFEMethod_DriverComplexForLazy.hpp"
#include "ElasticFEMethod.hpp"

#include "complex_for_lazy.h"

namespace MoFEM {

struct SetPositionsEntMethod: public moabField::EntMethod {
    ErrorCode rval;
    PetscErrorCode ierr;
    Interface& moab;

    EntityHandle node;
    double coords[3];

    SetPositionsEntMethod(Interface& _moab): EntMethod(),moab(_moab),node(no_handle) {}
    
    PetscErrorCode preProcess() {
      PetscFunctionBegin;
      PetscPrintf(PETSC_COMM_WORLD,"Start Set Positions\n");
      PetscFunctionReturn(0);
    } 
     
    PetscErrorCode operator()() {
      PetscFunctionBegin;
      if(dof_ptr->get_ent_type()!=MBVERTEX) PetscFunctionReturn(0);
      EntityHandle ent = dof_ptr->get_ent();
      int dof_rank = dof_ptr->get_dof_rank();
      double &fval = dof_ptr->get_FieldData();
      if(node!=ent) {
	rval = moab.get_coords(&ent,1,coords); CHKERR_PETSC(rval);
	node = ent;
      }
      fval = coords[dof_rank];
      PetscFunctionReturn(0);
    }

    PetscErrorCode postProcess() {
      PetscFunctionBegin;
      PetscPrintf(PETSC_COMM_WORLD,"End Set Positions\n");
      PetscFunctionReturn(0);
    }

};

struct NL_ElasticFEMethod: public FEMethod_DriverComplexForLazy {

  Range& SideSet2;

  NL_ElasticFEMethod(Interface& _moab,BaseDirihletBC *_dirihlet_bc_method_ptr,double _lambda,double _mu,Range &_SideSet2,int _verbose = 0): 
      FEMethod_DriverComplexForLazy(_moab,_dirihlet_bc_method_ptr,_lambda,_mu,_verbose), SideSet2(_SideSet2)  {

    set_PhysicalEquationNumber(hooke);
    //set_PhysicalEquationNumber(neohookean);

  }

  PetscErrorCode operator()() {
    PetscFunctionBegin;
    ierr = FEMethod_DriverComplexForLazy::operator()(SideSet2); CHKERRQ(ierr);
    PetscFunctionReturn(0);
  }


};

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
  C_SURFACE_FEMethod(Interface& _moab,EntityHandle skin_faces_meshset,Mat _C,int _verbose = 0): 
    FEMethod(),moab(_moab),C(_C) {
    diffNTRI.resize(6);
    ShapeDiffMBTRI(&diffNTRI[0]);
    g_NTRI3.resize(3*3);
    ShapeMBTRI(&g_NTRI3[0],G_TRI_X3,G_TRI_Y3,3);
    G_TRI_W = G_TRI_W3;
    double def_VAL[3*9];
    fill(&def_VAL[0],&def_VAL[3*9],0);
    rval = moab.tag_get_handle("MATERIAL_NORMAL",3,MB_TYPE_DOUBLE,th_material_normal,MB_TAG_CREAT|MB_TAG_SPARSE,&def_VAL); CHKERR_THROW(rval);
    rval = moab.get_entities_by_type(skin_faces_meshset,MBTRI,skin_faces,true);  CHKERR_THROW(rval);
  }

  PetscErrorCode preProcess() {
    PetscFunctionBegin;
    PetscFunctionReturn(0);
  }

  ublas::matrix<double> C_MAT_ELEM;
  ublas::vector<DofIdx> ent_global_col_indices,ent_global_row_indices;
  ublas::vector<double,ublas::bounded_array<double,9> > ent_dofs_data;
  ublas::vector<double,ublas::bounded_array<double,3> > ent_normal_map;
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
      ent_normal_map.resize(3);
      ierr = ShapeFaceNormalMBTRI(&diffNTRI[0],&ent_dofs_data.data()[0],&ent_normal_map.data()[0]); CHKERRQ(ierr);
      ent_normal_map *= siit->sense;
      //double area = norm_2(ent_normal_map);
      rval = moab.tag_set_data(th_material_normal,&face,1,&ent_normal_map.data()[0]); CHKERR_PETSC(rval);
      C_MAT_ELEM.resize(3,9);
      ublas::noalias(C_MAT_ELEM) = ublas::zero_matrix<double>(3,9);
      for(int gg = 0;gg<g_NTRI3.size()/3;gg++) {
	for(int nn = 0;nn<3;nn++) {
	  for(int dd = 0;dd<3;dd++) {
	    for(int nnn = 0;nnn<3;nnn++) {
	      C_MAT_ELEM(nn,3*nnn+dd) += G_TRI_W[gg]*ent_normal_map[dd]*g_NTRI3[3*gg+nn]*g_NTRI3[3*gg+nnn];
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
    }
    PetscFunctionReturn(0);
  }
  PetscErrorCode postProcess() {
    PetscFunctionBegin;
    PetscFunctionReturn(0);
  }

};

struct C_EDGE_FEMethod:public moabField::FEMethod {
  ErrorCode rval;
  PetscErrorCode ierr;
  Interface& moab;

  Mat C;
  Range edges,faces;
  vector<double> diffNTRI;
  vector<double> g_NTRI3;
  const double *G_TRI_W;
  C_EDGE_FEMethod(Interface& _moab,EntityHandle edges_meshset,Mat _C,int _verbose = 0): 
    FEMethod(),moab(_moab),C(_C) {
    diffNTRI.resize(6);
    ShapeDiffMBTRI(&diffNTRI[0]);
    g_NTRI3.resize(3*3);
    ShapeMBTRI(&g_NTRI3[0],G_TRI_X3,G_TRI_Y3,3);
    G_TRI_W = G_TRI_W3;
    rval = moab.get_entities_by_type(edges_meshset,MBEDGE,edges,true);  CHKERR_THROW(rval);
    rval = moab.get_entities_by_type(edges_meshset,MBTRI,faces,true);  CHKERR_THROW(rval);
  }

  map<EntityHandle,vector<EntityHandle> > edges_face_map;
  PetscErrorCode preProcess() {
    PetscFunctionBegin;
    edges_face_map.clear();
    PetscFunctionReturn(0);
  }

  ublas::matrix<double> C_MAT_ELEM;
  ublas::vector<DofIdx> ent_global_col_indices,ent_global_row_indices;
  ublas::vector<double,ublas::bounded_array<double,9> > ent_dofs_data;
  ublas::vector<double,ublas::bounded_array<double,3> > ent_normal_map;
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
	if(mit==edges_face_map.end()) {
	  edges_face_map[*eit].push_back(face);
	  eq_nb = 0;
	} else {
	  if(mit->second.size()==1) {
	    if(mit->second[0]==face) eq_nb = 0;
	    else mit->second.push_back(face);
	    eq_nb = 1;
	  } else {
	    if(mit->second[0]==face) eq_nb = 0;
	    else if(mit->second[1]==face) eq_nb = 1;
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
	    ent_global_row_indices[nn] = -1;
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
	ent_normal_map.resize(3);
	ierr = ShapeFaceNormalMBTRI(&diffNTRI[0],&ent_dofs_data.data()[0],&ent_normal_map.data()[0]); CHKERRQ(ierr);
	ent_normal_map *= siit->sense;
	C_MAT_ELEM.resize(3,9);
	ublas::noalias(C_MAT_ELEM) = ublas::zero_matrix<double>(3,9);
	for(int gg = 0;gg<g_NTRI3.size()/3;gg++) {
	  for(int nn = 0;nn<3;nn++) {
	    for(int dd = 0;dd<3;dd++) {
	      for(int nnn = 0;nnn<3;nnn++) {
		C_MAT_ELEM(nn,3*nnn+dd) += G_TRI_W[gg]*ent_normal_map[dd]*g_NTRI3[3*gg+nn]*g_NTRI3[3*gg+nnn];
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
      }
    }
    PetscFunctionReturn(0);
  }
  PetscErrorCode postProcess() {
    PetscFunctionBegin;
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
    }
    PetscFunctionReturn(0);
  }
  PetscErrorCode postProcess() {
    PetscFunctionBegin;
    PetscFunctionReturn(0);
  }

};


struct matPROJ_ctx {
  moabField& mField;
  Mat C,CT,CCT;
  KSP ksp;
  Vec _x_;
  VecScatter scatter;
  Vec Cx,CCTm1_Cx,CT_CCTm1_Cx;
  string x_problem,y_problem;
  bool init;
  bool debug;
  matPROJ_ctx(moabField& _mField,Mat _C,Mat _CT,Mat _CCT,string _x_problem,string _y_problem): mField(_mField),C(_C),CT(_CT),CCT(_CCT),
    x_problem(_x_problem),y_problem(_y_problem),init(true),debug(true) {}
  PetscErrorCode Destroy() {
    PetscFunctionBegin;
    if(init) PetscFunctionReturn(0);
    PetscErrorCode ierr;
    ierr = KSPDestroy(&ksp); CHKERRQ(ierr);
    ierr = VecScatterDestroy(&scatter); CHKERRQ(ierr);
    ierr = VecDestroy(&_x_); CHKERRQ(ierr);
    ierr = VecDestroy(&Cx); CHKERRQ(ierr);
    ierr = VecDestroy(&CCTm1_Cx); CHKERRQ(ierr);
    ierr = VecDestroy(&CT_CCTm1_Cx); CHKERRQ(ierr);
    PetscFunctionReturn(0);
  }
  friend PetscErrorCode matQ_mult_shell(Mat Q,Vec x,Vec f);
  friend PetscErrorCode matP_mult_shell(Mat P,Vec x,Vec f);
};

PetscErrorCode matQ_mult_shell(Mat Q,Vec x,Vec f) {
  PetscFunctionBegin;
  PetscErrorCode ierr;
  void *void_ctx;
  ierr = MatShellGetContext(Q,&void_ctx); CHKERRQ(ierr);
  matPROJ_ctx *ctx = (matPROJ_ctx*)void_ctx;
  if(ctx->init) {
    ctx->init = false;
    ierr = KSPCreate(PETSC_COMM_WORLD,&(ctx->ksp)); CHKERRQ(ierr);
    ierr = KSPSetOperators(ctx->ksp,ctx->CCT,ctx->CCT,SAME_NONZERO_PATTERN); CHKERRQ(ierr);
    ierr = KSPSetFromOptions(ctx->ksp); CHKERRQ(ierr);
    ierr = KSPSetUp(ctx->ksp); CHKERRQ(ierr);
    ierr = MatGetVecs(ctx->C,&ctx->_x_,PETSC_NULL); CHKERRQ(ierr);
    ierr = MatGetVecs(ctx->C,PETSC_NULL,&ctx->Cx); CHKERRQ(ierr);
    ierr = MatGetVecs(ctx->CCT,PETSC_NULL,&ctx->CCTm1_Cx); CHKERRQ(ierr);
    ierr = VecDuplicate(ctx->_x_,&ctx->CT_CCTm1_Cx); CHKERRQ(ierr);
    ierr = ctx->mField.VecScatterCreate(x,ctx->x_problem,Row,ctx->_x_,ctx->y_problem,Col,&ctx->scatter); CHKERRQ(ierr);
  }
  ierr = VecCopy(x,f); CHKERRQ(ierr);
  ierr = VecScatterBegin(ctx->scatter,x,ctx->_x_,INSERT_VALUES,SCATTER_FORWARD); CHKERRQ(ierr);
  ierr = VecScatterEnd(ctx->scatter,x,ctx->_x_,INSERT_VALUES,SCATTER_FORWARD); CHKERRQ(ierr);
  //ierr = VecView(ctx->_x_,PETSC_VIEWER_STDOUT_WORLD); CHKERRQ(ierr);
  if(ctx->debug) {
    ierr = VecScatterBegin(ctx->scatter,ctx->_x_,f,INSERT_VALUES,SCATTER_REVERSE); CHKERRQ(ierr);
    ierr = VecScatterEnd(ctx->scatter,ctx->_x_,f,INSERT_VALUES,SCATTER_REVERSE); CHKERRQ(ierr);
    PetscBool  flg;
    ierr = VecEqual(x,f,&flg); CHKERRQ(ierr);
    if(flg ==  PETSC_FALSE) SETERRQ(PETSC_COMM_SELF,1,"scatter is not working");
  }
  ierr = MatMult(ctx->C,ctx->_x_,ctx->Cx);  CHKERRQ(ierr);
  ierr = KSPSolve(ctx->ksp,ctx->Cx,ctx->CCTm1_Cx); CHKERRQ(ierr);
  ierr = MatMult(ctx->CT,ctx->CCTm1_Cx,ctx->CT_CCTm1_Cx);  CHKERRQ(ierr);
  ierr = VecScale(ctx->CT_CCTm1_Cx,-1); CHKERRQ(ierr);
  ierr = VecScatterBegin(ctx->scatter,ctx->CT_CCTm1_Cx,f,ADD_VALUES,SCATTER_REVERSE); CHKERRQ(ierr);
  ierr = VecScatterEnd(ctx->scatter,ctx->CT_CCTm1_Cx,f,ADD_VALUES,SCATTER_REVERSE); CHKERRQ(ierr);
  PetscFunctionReturn(0);
}
PetscErrorCode matP_mult_shell(Mat P,Vec x,Vec f) {
  PetscFunctionBegin;
  PetscErrorCode ierr;
  void *void_ctx;
  ierr = MatShellGetContext(P,&void_ctx); CHKERRQ(ierr);
  matPROJ_ctx *ctx = (matPROJ_ctx*)void_ctx;
  if(ctx->init) {
    ctx->init = false;
    ierr = KSPCreate(PETSC_COMM_WORLD,&(ctx->ksp)); CHKERRQ(ierr);
    ierr = KSPSetOperators(ctx->ksp,ctx->CCT,ctx->CCT,SAME_NONZERO_PATTERN); CHKERRQ(ierr);
    ierr = KSPSetFromOptions(ctx->ksp); CHKERRQ(ierr);
    ierr = KSPSetUp(ctx->ksp); CHKERRQ(ierr);
    ierr = MatGetVecs(ctx->C,&ctx->_x_,PETSC_NULL); CHKERRQ(ierr);
    ierr = MatGetVecs(ctx->C,PETSC_NULL,&ctx->Cx); CHKERRQ(ierr);
    ierr = MatGetVecs(ctx->CCT,PETSC_NULL,&ctx->CCTm1_Cx); CHKERRQ(ierr);
    ierr = VecDuplicate(ctx->_x_,&ctx->CT_CCTm1_Cx); CHKERRQ(ierr);
    ierr = ctx->mField.VecScatterCreate(x,ctx->x_problem,Row,ctx->_x_,ctx->y_problem,Col,&ctx->scatter); CHKERRQ(ierr);
  }
  ierr = VecScatterBegin(ctx->scatter,x,ctx->_x_,INSERT_VALUES,SCATTER_FORWARD); CHKERRQ(ierr);
  ierr = VecScatterEnd(ctx->scatter,x,ctx->_x_,INSERT_VALUES,SCATTER_FORWARD); CHKERRQ(ierr);
  ierr = MatMult(ctx->C,ctx->_x_,ctx->Cx);  CHKERRQ(ierr);
  ierr = KSPSolve(ctx->ksp,ctx->Cx,ctx->CCTm1_Cx); CHKERRQ(ierr);
  ierr = MatMult(ctx->CT,ctx->CCTm1_Cx,ctx->CT_CCTm1_Cx);  CHKERRQ(ierr);
  ierr = VecZeroEntries(f); CHKERRQ(ierr);
  ierr = VecGhostUpdateBegin(f,INSERT_VALUES,SCATTER_FORWARD); CHKERRQ(ierr);
  ierr = VecGhostUpdateEnd(f,INSERT_VALUES,SCATTER_FORWARD); CHKERRQ(ierr);
  ierr = VecScatterBegin(ctx->scatter,ctx->CT_CCTm1_Cx,f,ADD_VALUES,SCATTER_REVERSE); CHKERRQ(ierr);
  ierr = VecScatterEnd(ctx->scatter,ctx->CT_CCTm1_Cx,f,ADD_VALUES,SCATTER_REVERSE); CHKERRQ(ierr);
  PetscFunctionReturn(0);
}



struct MaterialForcesFEMethod: public FEMethod_DriverComplexForLazy {

  Vec F_MATERIAL;
  MaterialForcesFEMethod(Interface& _moab,BaseDirihletBC *_dirihlet_bc_method_ptr,double _lambda,double _mu,Vec _F_MATERIAL,int _verbose = 0): 
      FEMethod_DriverComplexForLazy(_moab,_dirihlet_bc_method_ptr,_lambda,_mu,_verbose),F_MATERIAL(_F_MATERIAL) {

    set_PhysicalEquationNumber(hooke);
    type_of_analysis = material_analysis;

  }

  PetscErrorCode preProcess() {
    PetscFunctionBegin;
    PetscFunctionReturn(0);
  }

  PetscErrorCode operator()() {
    PetscFunctionBegin;
    ierr = OpComplexForLazyStart(); CHKERRQ(ierr);
    ierr = GetIndices(RowGlobMaterial,ColGlobMaterial,material_field_name); CHKERRQ(ierr);
    ierr = GetData(dofs_x_edge_data,dofs_x_edge,
      dofs_x_face_data,dofs_x_face,
      dofs_x_volume,dofs_x,
      spatial_field_name); CHKERRQ(ierr);
    ierr = GetFint(); CHKERRQ(ierr);
    ierr = VecSetValues(F_MATERIAL,RowGlobMaterial[0].size(),&(RowGlobMaterial[0])[0],&(Fint_H.data())[0],ADD_VALUES); CHKERRQ(ierr);
    PetscFunctionReturn(0);
  }

  PetscErrorCode postProcess() {
    PetscFunctionBegin;
    PetscFunctionReturn(0);
  }


};

}

#endif //__NONLINEAR_ELASTICITY_HPP__
