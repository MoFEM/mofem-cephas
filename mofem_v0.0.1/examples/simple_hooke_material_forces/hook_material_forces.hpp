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

struct C_MATRIX_FEMethod:public moabField::FEMethod {
  ErrorCode rval;
  PetscErrorCode ierr;
  Interface& moab;

  Mat C;
  Range skin_faces;
  Tag th_material_normal;
  Tag th_c_mat_elem[4];
  vector<double> diffNTRI;
  vector<double> g_NTRI3;
  const double *G_TRI_W;
  C_MATRIX_FEMethod(Interface& _moab,EntityHandle skin_faces_meshset,Mat _C,int _verbose = 0): FEMethod(),moab(_moab),C(_C) {
    diffNTRI.resize(6);
    ShapeDiffMBTRI(&diffNTRI[0]);
    g_NTRI3.resize(3*3);
    ShapeMBTRI(&g_NTRI3[0],G_TRI_X3,G_TRI_Y3,3);
    G_TRI_W = G_TRI_W3;
    double def_VAL[3*9];
    fill(&def_VAL[0],&def_VAL[3*9],0);
    rval = moab.tag_get_handle("MATERIAL_NORMAL",3,MB_TYPE_DOUBLE,th_material_normal,MB_TAG_CREAT|MB_TAG_SPARSE,&def_VAL); CHKERR_THROW(rval);
    rval = moab.tag_get_handle("C_MAT_ELEM_SIDE0",3*9,MB_TYPE_DOUBLE,th_c_mat_elem[0],MB_TAG_CREAT|MB_TAG_SPARSE,&def_VAL); CHKERR_THROW(rval);
    rval = moab.tag_get_handle("C_MAT_ELEM_SIDE1",3*9,MB_TYPE_DOUBLE,th_c_mat_elem[1],MB_TAG_CREAT|MB_TAG_SPARSE,&def_VAL); CHKERR_THROW(rval);
    rval = moab.tag_get_handle("C_MAT_ELEM_SIDE2",3*9,MB_TYPE_DOUBLE,th_c_mat_elem[2],MB_TAG_CREAT|MB_TAG_SPARSE,&def_VAL); CHKERR_THROW(rval);
    rval = moab.tag_get_handle("C_MAT_ELEM_SIDE3",3*9,MB_TYPE_DOUBLE,th_c_mat_elem[3],MB_TAG_CREAT|MB_TAG_SPARSE,&def_VAL); CHKERR_THROW(rval);
    rval = moab.get_entities_by_type(skin_faces_meshset,MBTRI,skin_faces,true);  CHKERR_THROW(rval);
  }

  PetscErrorCode preProcess() {
    PetscFunctionBegin;
    PetscFunctionReturn(0);
  }

  ublas::vector<double> C_CONST_SHAPE;
  ublas::matrix<double> C_MAT_ELEM;
  vector<DofIdx> ent_global_col_indices,ent_global_row_indices;
  ublas::vector<double,ublas::bounded_array<double,9> > ent_dofs_data;
  ublas::vector<double,ublas::bounded_array<double,3> > ent_normal_map;
  PetscErrorCode operator()() {
    PetscFunctionBegin;
    SideNumber_multiIndex &side_table = fe_ptr->get_side_number_table();
    SideNumber_multiIndex::nth_index<1>::type::iterator siit = side_table.get<1>().lower_bound(boost::make_tuple(MBTRI,0));
    SideNumber_multiIndex::nth_index<1>::type::iterator hi_siit = side_table.get<1>().upper_bound(boost::make_tuple(MBTRI,4));
    for(;siit!=hi_siit;siit++) {
      EntityHandle face = siit->ent;
      if(skin_faces.end()==find(skin_faces.begin(),skin_faces.end(),face)) continue;
      ent_dofs_data.resize(9);
      ent_global_col_indices.resize(9);
      ent_global_row_indices.resize(3);
      const EntityHandle* conn_face; 
      int num_nodes; 
      rval = moab.get_connectivity(face,conn_face,num_nodes,true); CHKERR_PETSC(rval);
      for(int nn = 0;nn<num_nodes;nn++) {
	FENumeredDofMoFEMEntity_multiIndex::index<Composite_mi_tag3>::type::iterator dit,hi_dit;
	dit = col_multiIndex->get<Composite_mi_tag3>().lower_bound(boost::make_tuple("MESH_NODE_POSITIONS",conn_face[nn]));
	hi_dit = col_multiIndex->get<Composite_mi_tag3>().upper_bound(boost::make_tuple("MESH_NODE_POSITIONS",conn_face[nn]));
	if(distance(dit,hi_dit)!=3) SETERRQ(PETSC_COMM_SELF,1,"data inconsistency, number of dof on node for MESH_NODE_POSITIONS should be 3");
	for(;dit!=hi_dit;dit++) {
	  ent_dofs_data[nn*3+dit->get_dof_rank()] = dit->get_FieldData();
	  int global_idx = dit->get_petsc_gloabl_dof_idx();
	  ent_global_col_indices[nn*3+dit->get_dof_rank()] = global_idx;
	  int local_idx = dit->get_petsc_local_dof_idx();
	  if(local_idx<0) SETERRQ(PETSC_COMM_SELF,1,"data inconsistency");
	}
	dit = row_multiIndex->get<Composite_mi_tag3>().lower_bound(boost::make_tuple("CONST_SHAPE_LAMBDA",conn_face[nn]));
	hi_dit = row_multiIndex->get<Composite_mi_tag3>().upper_bound(boost::make_tuple("CONST_SHAPE_LAMBDA",conn_face[nn]));
	if(distance(dit,hi_dit)!=1) SETERRQ(PETSC_COMM_SELF,1,"data inconsistency, number of dof on node for CONST_SHAPE_LAMBDA should be 1");
	int global_idx = dit->get_petsc_gloabl_dof_idx();
	ent_global_row_indices[nn] = global_idx;
	int local_idx = dit->get_petsc_local_dof_idx();
	if(local_idx<0) SETERRQ(PETSC_COMM_SELF,1,"data inconsistency");
      }
      ent_normal_map.resize(3);
      ierr = ShapeFaceNormalMBTRI(&diffNTRI[0],&ent_dofs_data.data()[0],&ent_normal_map.data()[0]); CHKERRQ(ierr);
      rval = moab.tag_set_data(th_material_normal,&face,1,&ent_normal_map.data()[0]); CHKERR_PETSC(rval);
      C_CONST_SHAPE.resize(3*3);
      for(int nn = 0;nn<num_nodes;nn++) {
	for(int dd = 0;dd<3;dd++) {
	  C_CONST_SHAPE[nn*3+dd] = ent_normal_map[dd];
	}
      }
      //C_CONST_SHAPE /= norm_2(ent_normal_map);
      C_MAT_ELEM.resize(3,9);
      for(int gg = 0;gg<g_NTRI3.size()/3;gg++) {
	for(int nn = 0;nn<3;nn++) {
	  for(int dd = 0;dd<9;dd++) {
	    if(gg == 0) C_MAT_ELEM(nn,dd) = G_TRI_W[gg]*C_CONST_SHAPE[dd]*g_NTRI3[3*gg+nn];
	    else C_MAT_ELEM(nn,dd) += G_TRI_W[gg]*C_CONST_SHAPE[dd]*g_NTRI3[3*gg+nn];
	  }
	}
      }
      EntityHandle fe_ent = fe_ptr->get_ent();
      int side_number = siit->side_number;
      rval = moab.tag_set_data(th_c_mat_elem[side_number],&fe_ent,1,&(C_MAT_ELEM.data()[0])); CHKERR_PETSC(rval);
      ierr = MatSetValues(C,
	ent_global_row_indices.size(),&(ent_global_row_indices[0]),
	ent_global_col_indices.size(),&(ent_global_col_indices[0]),
	&(C_MAT_ELEM.data())[0],ADD_VALUES); CHKERRQ(ierr);
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
  matPROJ_ctx(moabField& _mField,Mat _C,Mat _CT,Mat _CCT): mField(_mField),C(_C),CT(_CT),CCT(_CCT),
    x_problem("MATERIAL_MECHANICS"),y_problem("C_MATRIX"),init(true) {}
  friend PetscErrorCode matQTAQ_mult_shell(Mat QTAQ,Vec x,Vec f);
};

PetscErrorCode matQTAQ_mult_shell(Mat QTAQ,Vec x,Vec f) {
  PetscFunctionBegin;
  PetscErrorCode ierr;
  void *void_ctx;
  ierr = MatShellGetContext(QTAQ,&void_ctx); CHKERRQ(ierr);
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
    ierr = ctx->mField.VecScatterCreate(x,ctx->x_problem,Col,ctx->_x_,ctx->y_problem,Col,&ctx->scatter,2); CHKERRQ(ierr);
  }
  ierr = VecScatterBegin(ctx->scatter,x,ctx->_x_,INSERT_VALUES,SCATTER_FORWARD); CHKERRQ(ierr);
  ierr = VecScatterEnd(ctx->scatter,x,ctx->_x_,INSERT_VALUES,SCATTER_FORWARD); CHKERRQ(ierr);
  ierr = MatMult(ctx->C,ctx->_x_,ctx->Cx);  CHKERRQ(ierr);
  ierr = KSPSolve(ctx->ksp,ctx->Cx,ctx->CCTm1_Cx); CHKERRQ(ierr);
  ierr = MatMult(ctx->CT,ctx->CCTm1_Cx,ctx->CT_CCTm1_Cx);  CHKERRQ(ierr);
  ierr = VecScatterBegin(ctx->scatter,f,ctx->CT_CCTm1_Cx,INSERT_VALUES,SCATTER_REVERSE); CHKERRQ(ierr);
  ierr = VecScatterEnd(ctx->scatter,f,ctx->CT_CCTm1_Cx,INSERT_VALUES,SCATTER_REVERSE); CHKERRQ(ierr);
  ierr = VecAYPX(f,-1,x); CHKERRQ(ierr);
  ierr = VecScale(f,-1); CHKERRQ(ierr);
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
