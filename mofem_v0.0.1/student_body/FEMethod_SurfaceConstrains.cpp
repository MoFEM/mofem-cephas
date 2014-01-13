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

#include "FEMethod_SurfaceConstrains.hpp"

#include <math.h>
#include <complex>
extern "C" {
#include <complex.h>
void tetcircumcenter_tp(double a[3],double b[3],double c[3], double d[3],
  double circumcenter[3],double *xi,double *eta,double *zeta);
void tricircumcenter3d_tp(double a[3],double b[3],double c[3],
  double circumcenter[3],double *xi,double *eta);
}

namespace MoFEM {

const int debug_constrains = 1;

PetscErrorCode C_SURFACE_FEMethod::preProcess() {
    PetscFunctionBegin;
    PetscFunctionReturn(0);
  }

void C_SURFACE_FEMethod::run_in_constructor() {
    diffNTRI.resize(6);
    ShapeDiffMBTRI(&diffNTRI[0]);
    g_NTRI.resize(7*3);
    ShapeMBTRI(&g_NTRI[0],G_TRI_X7,G_TRI_Y7,7);
    G_TRI_W = G_TRI_W7;
    double def_VAL[3*9];
    fill(&def_VAL[0],&def_VAL[3*9],0);
    lambda_global_row_indices.resize(3);
    lambda_global_col_indices.resize(3);
    dof_global_row_indices.resize(1+6+4+1);
    dof_global_col_indices.resize(1+6+4+1);
    dof_global_row_indices[0].resize(9);
    dof_global_col_indices[0].resize(9);
    ent_lambda_data.resize(3);
    ent_dofs_data.resize(9);
    ent_idofs_data.resize(9);
    coords.resize(9);
    //
    C_MAT_ELEM.resize(3,9);
    iC_MAT_ELEM.resize(3,9);
    CT_MAT_ELEM.resize(9,3);
    dC_MAT_ELEM.resize(3,9);
    dCT_MAT_ELEM.resize(9,9);
  }
 
C_SURFACE_FEMethod::C_SURFACE_FEMethod(Interface& _moab,BaseDirihletBC *_dirihlet_bc_method_ptr,Mat _C,string _lambda_field_name,int _verbose): 
    FEMethod(),moab(_moab),dirihlet_bc_method_ptr(_dirihlet_bc_method_ptr),C(_C),lambda_field_name(_lambda_field_name),updated(false) {
    run_in_constructor();
  }
C_SURFACE_FEMethod::C_SURFACE_FEMethod(Interface& _moab,BaseDirihletBC *_dirihlet_bc_method_ptr,Mat _C,int _verbose): 
    FEMethod(),moab(_moab),dirihlet_bc_method_ptr(_dirihlet_bc_method_ptr),C(_C),lambda_field_name("LAMBDA_SURFACE"),updated(false) {
    run_in_constructor();
  }

PetscErrorCode C_SURFACE_FEMethod::cOnstrain(double *dofs_X,double *dofs_iX,double *C,double *iC) {
  PetscFunctionBegin;
  //set complex material position vector
  __CLPK_doublecomplex x_dofs_X[9];
  for(int nn = 0;nn<3;nn++) {
    for(int dd = 0;dd<3;dd++) {
      x_dofs_X[nn*3+dd].r = dofs_X[nn*3+dd];
      if(dofs_iX != NULL) {
	x_dofs_X[nn*3+dd].i = dofs_iX[nn*3+dd];
      } else {
	x_dofs_X[nn*3+dd].i = 0;
      }
  }}
  //calulate normal
  __CLPK_doublecomplex x_normal[3];
  ierr = ShapeFaceNormalMBTRI_complex(&diffNTRI[0],x_dofs_X,x_normal); CHKERRQ(ierr);
  //set direction if crack or interface surface
  Tag th_internal_node;
  const EntityHandle def_node[] = {0};
  rval = moab.tag_get_handle("INTERNAL_SIDE_NODE",1,MB_TYPE_HANDLE,
      th_internal_node,MB_TAG_CREAT|MB_TAG_SPARSE,def_node);
  EntityHandle internal_node;
  rval = moab.tag_get_data(th_internal_node,&face,1,&internal_node); CHKERR_PETSC(rval);
  if(internal_node!=0) {
    Tag th_interface_side;
    rval = moab.tag_get_handle("INTERFACE_SIDE",th_interface_side); CHKERR_PETSC(rval);
    int side;
    rval = moab.tag_get_data(th_interface_side,&face,1,&side); CHKERR_PETSC(rval);
    double coords_internal_node[3];
    rval = moab.get_coords(&internal_node,1,coords_internal_node); CHKERR_PETSC(rval);
    cblas_daxpy(3,-1,dofs_X,1,coords_internal_node,1);
    if(side) cblas_dscal(3,-1,coords_internal_node,1);
    __CLPK_doublecomplex xdot = { 
      copysign(1,x_normal[0].r*coords_internal_node[0]+
      x_normal[1].r*coords_internal_node[1]+
      x_normal[2].r*coords_internal_node[2]), 0 };
    cblas_zscal(3,&xdot,x_normal,1);
  };
  //calulare complex normal length
  /*double __complex__ xarea = csqrt(
      cpow((x_normal[0].r+I*x_normal[0].i),2)+
      cpow((x_normal[1].r+I*x_normal[1].i),2)+
      cpow((x_normal[2].r+I*x_normal[2].i),2));
  double normal[3] = { x_normal[0].r/creal(xarea), x_normal[1].r/creal(xarea), x_normal[2].r/creal(xarea) };*/
  /*const int def_normal[3] = { 0,0,0 };
  Tag th_normal;
  rval = moab.tag_get_handle("NORMAL_TEST",3,MB_TYPE_DOUBLE,
    th_normal,MB_TAG_CREAT|MB_TAG_SPARSE,def_normal);
  rval = moab.tag_set_data(th_normal,&face,1,normal); CHKERR_PETSC(rval);*/
  //scale normal vector
  for(int dd = 0;dd<3;dd++) {
    double __complex__ val;
    val = (x_normal[dd].r+I*x_normal[dd].i);///xarea;
    x_normal[dd].r = creal(val);
    x_normal[dd].i = cimag(val);
  }
  if( C!=NULL) bzero( C,3*9*sizeof(double));
  if(iC!=NULL) bzero(iC,3*9*sizeof(double));
  for(unsigned int gg = 0;gg<g_NTRI.size()/3;gg++) {
    //->row
    for(int nn = 0;nn<3;nn++) {
      if(lambda_global_row_indices[nn]==-1) continue;
      //->col
      for(int mm = 0;mm<3;mm++) {
	if(lambda_global_row_indices[mm]==-1) continue;
	for(int dd = 0;dd<3;dd++) {
	  if( C!=NULL) {
	     C[9*nn + 3*mm+dd] += G_TRI_W[gg]*( g_NTRI[3*gg+nn]*g_NTRI[3*gg+mm]*x_normal[dd].r );
	  }
	  if(iC!=NULL) {
	    iC[9*nn + 3*mm+dd] += G_TRI_W[gg]*( g_NTRI[3*gg+nn]*g_NTRI[3*gg+mm]*x_normal[dd].i );
	  }
	}
      }
      //<-col
    }
    //<-row
  }
  PetscFunctionReturn(0);
}
PetscErrorCode C_SURFACE_FEMethod::iNtegrate(bool transpose,bool nonlinear) {
  PetscFunctionBegin;
  const double eps = 1e-10;
  double center[3]; 
  tricircumcenter3d_tp(&coords.data()[0],&coords.data()[3],&coords.data()[6],center,NULL,NULL);
  cblas_daxpy(3,-1,&coords.data()[0],1,center,1);
  double r = cblas_dnrm2(3,center,1);
  try {
    if(nonlinear||updated) {
      ierr = cOnstrain(&*ent_dofs_data.data().begin(),NULL,&*C_MAT_ELEM.data().begin(),NULL); CHKERRQ(ierr);
    } else {
      ierr = cOnstrain(&*coords.data().begin(),NULL,&*C_MAT_ELEM.data().begin(),NULL); CHKERRQ(ierr);
    }
    if(transpose) {
      ublas::noalias(CT_MAT_ELEM) = trans(C_MAT_ELEM);
      if(nonlinear) {
	ublas::noalias(dC_MAT_ELEM) = ublas::zero_matrix<double>(3,9);
	ublas::noalias(dCT_MAT_ELEM) = ublas::zero_matrix<double>(9,9);
	for(int dd = 0;dd<9;dd++) {
	  ublas::noalias(ent_idofs_data) = ublas::zero_vector<double>(9);
	  ent_idofs_data[dd] = r*eps;
	  ierr = cOnstrain(&*ent_dofs_data.data().begin(),&*ent_idofs_data.data().begin(),NULL,&*iC_MAT_ELEM.data().begin()); CHKERRQ(ierr);
	  ig_VEC_ELEM = prod(iC_MAT_ELEM/(r*eps),ent_dofs_data-coords);
	  for(int nnn = 0;nnn<3;nnn++) {
	    dC_MAT_ELEM(nnn,dd) += ig_VEC_ELEM[nnn];
	  }
	  if(transpose) {
	    if_VEC_ELEM = prod(trans(iC_MAT_ELEM)/(r*eps),ent_lambda_data);
	    for(int nnn = 0;nnn<3;nnn++) {
	      for(int ddd = 0;ddd<3;ddd++) {
		dCT_MAT_ELEM(3*nnn+ddd,dd) += if_VEC_ELEM[3*nnn+ddd];
	      }
	    }
	  }
	}
      }
    }
  } catch (const std::exception& ex) {
    ostringstream ss;
    ss << "thorw in method: " << ex.what() << " at line " << __LINE__ << " in file " << __FILE__ << endl;
    SETERRQ(PETSC_COMM_SELF,1,ss.str().c_str());
  }
  PetscFunctionReturn(0);
}

PetscErrorCode C_SURFACE_FEMethod::aSsemble(bool transpose,bool nonlinear) {
    PetscFunctionBegin;
    ierr = MatSetValues(C,
      lambda_global_row_indices.size(),&(lambda_global_row_indices.data()[0]),
      dof_global_col_indices[0].size(),&(dof_global_col_indices[0].data()[0]),
      &(C_MAT_ELEM.data())[0],ADD_VALUES); CHKERRQ(ierr);
    if(nonlinear) {
      ierr = MatSetValues(C,
	lambda_global_row_indices.size(),&(lambda_global_row_indices.data()[0]),
	dof_global_col_indices[0].size(),&(dof_global_col_indices[0].data()[0]),
	&(dC_MAT_ELEM.data())[0],ADD_VALUES); CHKERRQ(ierr);
    }
    if(transpose) {
      ierr = MatSetValues(C,
	dof_global_row_indices[0].size(),&(dof_global_row_indices[0].data()[0]),
	lambda_global_col_indices.size(),&(lambda_global_col_indices.data()[0]),
	&CT_MAT_ELEM.data()[0],ADD_VALUES); CHKERRQ(ierr);
      if(nonlinear) {
	ierr = MatSetValues(C,
	  dof_global_row_indices[0].size(),&(dof_global_row_indices[0].data()[0]),
	  dof_global_col_indices[0].size(),&(dof_global_col_indices[0].data()[0]),
	  &dCT_MAT_ELEM.data()[0],ADD_VALUES); CHKERRQ(ierr); 
      }
    }
    PetscFunctionReturn(0);
}

PetscErrorCode C_SURFACE_FEMethod::operator()(bool transpose,bool nonlinear) {
    PetscFunctionBegin;
    try {
      face = fe_ptr->get_ent();
      const EntityHandle* conn_face; 
      int num_nodes; 
      rval = moab.get_connectivity(face,conn_face,num_nodes,true); CHKERR_PETSC(rval);
      for(int nn = 0;nn<num_nodes;nn++) {
  	FENumeredDofMoFEMEntity_multiIndex::index<Composite_Name_And_Ent_mi_tag>::type::iterator dit,hi_dit;
  	dit = row_multiIndex->get<Composite_Name_And_Ent_mi_tag>().lower_bound(boost::make_tuple(lambda_field_name,conn_face[nn]));
  	hi_dit = row_multiIndex->get<Composite_Name_And_Ent_mi_tag>().upper_bound(boost::make_tuple(lambda_field_name,conn_face[nn]));
  	if(distance(dit,hi_dit)==0) {
  	  lambda_global_row_indices[nn] = -1;
	  ent_lambda_data[nn] = 0;
  	} else {
  	  if(distance(dit,hi_dit)!=1) SETERRQ1(PETSC_COMM_SELF,1,"data inconsistency, number of dof on node for < %s > should be 1",lambda_field_name.c_str());
	  ent_lambda_data[nn] = dit->get_FieldData();
  	  int global_idx = dit->get_petsc_gloabl_dof_idx();
  	  lambda_global_row_indices[nn] = global_idx;
  	  int local_idx = dit->get_petsc_local_dof_idx();
  	  if(local_idx<0) SETERRQ(PETSC_COMM_SELF,1,"data inconsistency, negative index of local dofs on element");
  	}
  	dit = col_multiIndex->get<Composite_Name_And_Ent_mi_tag>().lower_bound(boost::make_tuple("MESH_NODE_POSITIONS",conn_face[nn]));
  	hi_dit = col_multiIndex->get<Composite_Name_And_Ent_mi_tag>().upper_bound(boost::make_tuple("MESH_NODE_POSITIONS",conn_face[nn]));
  	if(distance(dit,hi_dit)!=3) SETERRQ(PETSC_COMM_SELF,1,"data inconsistency, number of dof on node for MESH_NODE_POSITIONS should be 3");
  	for(;dit!=hi_dit;dit++) {
  	  ent_dofs_data[nn*3+dit->get_dof_rank()] = dit->get_FieldData();
  	  int global_idx = dit->get_petsc_gloabl_dof_idx();
  	  assert(nn*3+dit->get_dof_rank()<9);
  	  dof_global_col_indices[0][nn*3+dit->get_dof_rank()] = global_idx;
  	  int local_idx = dit->get_petsc_local_dof_idx();
  	  if(local_idx<0) SETERRQ(PETSC_COMM_SELF,1,"data inconsistency, negative index of local dofs on element");
  	}
	if(transpose) {
	  dit = col_multiIndex->get<Composite_Name_And_Ent_mi_tag>().lower_bound(boost::make_tuple(lambda_field_name,conn_face[nn]));
	  hi_dit = col_multiIndex->get<Composite_Name_And_Ent_mi_tag>().upper_bound(boost::make_tuple(lambda_field_name,conn_face[nn]));
	  if(distance(dit,hi_dit)==0) {
	    lambda_global_col_indices[nn] = -1;
	  } else {
	    if(distance(dit,hi_dit)!=1) SETERRQ1(PETSC_COMM_SELF,1,"data inconsistency, number of dof on node for < %s > should be 1",lambda_field_name.c_str());
	    int global_idx = dit->get_petsc_gloabl_dof_idx();
	    lambda_global_col_indices[nn] = global_idx;
	    int local_idx = dit->get_petsc_local_dof_idx();
	    if(local_idx<0) {
	      SETERRQ(PETSC_COMM_SELF,1,"data inconsistency, negative index of local dofs on element");
	    }
	  }
	  dit = row_multiIndex->get<Composite_Name_And_Ent_mi_tag>().lower_bound(boost::make_tuple("MESH_NODE_POSITIONS",conn_face[nn]));
	  hi_dit = row_multiIndex->get<Composite_Name_And_Ent_mi_tag>().upper_bound(boost::make_tuple("MESH_NODE_POSITIONS",conn_face[nn]));
	  if(distance(dit,hi_dit)!=3) SETERRQ(PETSC_COMM_SELF,1,"data inconsistency, number of dof on node for MESH_NODE_POSITIONS should be 3");
	  for(;dit!=hi_dit;dit++) {
	    int global_idx = dit->get_petsc_gloabl_dof_idx();
	    assert(nn*3+dit->get_dof_rank()<9);
	    dof_global_row_indices[0][nn*3+dit->get_dof_rank()] = global_idx;
	    int local_idx = dit->get_petsc_local_dof_idx();
	    if(local_idx<0) {
	      SETERRQ(PETSC_COMM_SELF,1,"data inconsistency, negative index of local dofs on element");
	    }
	  }
	}
	DirihletBC.clear();
	ierr = dirihlet_bc_method_ptr->SetDirihletBC_to_ElementIndicies(
	  this,dof_global_row_indices,dof_global_col_indices,DirihletBC); CHKERRQ(ierr);
      }
      rval = moab.get_coords(conn_face,num_nodes,&*coords.data().begin()); CHKERR_PETSC(rval);
      ierr = this->iNtegrate(transpose,nonlinear); CHKERRQ(ierr);
      ierr = this->aSsemble(transpose,nonlinear); CHKERRQ(ierr);
    } catch (const std::exception& ex) {
      ostringstream ss;
      ss << "thorw in method: " << ex.what() << " at line " << __LINE__ << " in file " << __FILE__ << endl;
      SETERRQ(PETSC_COMM_SELF,1,ss.str().c_str());
    }
    PetscFunctionReturn(0);
  }

PetscErrorCode C_SURFACE_FEMethod::postProcess() {
  PetscFunctionBegin;
  PetscFunctionReturn(0);
}

g_SURFACE_FEMethod::g_SURFACE_FEMethod(Interface& _moab,BaseDirihletBC *_dirihlet_bc_method_ptr,Vec _g,string _lambda_field_name,int _verbose): 
    C_SURFACE_FEMethod(_moab,_dirihlet_bc_method_ptr,PETSC_NULL,_lambda_field_name,_verbose),g(_g) {
    g_VEC_ELEM.resize(3);
  }
g_SURFACE_FEMethod::g_SURFACE_FEMethod(Interface& _moab,BaseDirihletBC *_dirihlet_bc_method_ptr,Vec _g,int _verbose): 
    C_SURFACE_FEMethod(_moab,_dirihlet_bc_method_ptr,PETSC_NULL,_verbose),g(_g) {
    g_VEC_ELEM.resize(3);
  }

PetscErrorCode g_SURFACE_FEMethod::iNtegrate(bool transpose,bool nonlinear) {
  PetscFunctionBegin;
  try {
    C_SURFACE_FEMethod::iNtegrate(transpose,nonlinear);
    //
    g_VEC_ELEM = prod(C_MAT_ELEM,ent_dofs_data-coords);
    f_VEC_ELEM = prod(CT_MAT_ELEM,ent_lambda_data);
    //
    //ierr = cOnstrain(&*coords.data().begin(),NULL,&*C_MAT_ELEM.data().begin(),NULL); CHKERRQ(ierr);
    //g_VEC_ELEM = g_VEC_ELEM - prod(C_MAT_ELEM,coords);
  } catch (const std::exception& ex) {
    ostringstream ss;
    ss << "thorw in method: " << ex.what() << " at line " << __LINE__ << " in file " << __FILE__ << endl;
    SETERRQ(PETSC_COMM_SELF,1,ss.str().c_str());
  }
  PetscFunctionReturn(0);
}
PetscErrorCode g_SURFACE_FEMethod::aSsemble(bool transpose,bool nonlinear) {
  PetscFunctionBegin;
  if(lambda_global_row_indices.size()!=g_VEC_ELEM.size()) {
    SETERRQ2(PETSC_COMM_SELF,1,"data inconsistency %d != %d",
	lambda_global_row_indices.size(),g_VEC_ELEM.size());
  }
  ierr = VecSetOption(g,VEC_IGNORE_NEGATIVE_INDICES,PETSC_TRUE); CHKERRQ(ierr);
  ierr = VecSetValues(g,
    lambda_global_row_indices.size(),&(lambda_global_row_indices.data()[0]),
    &(g_VEC_ELEM.data())[0],ADD_VALUES); CHKERRQ(ierr);
  if(transpose) {
    if(dof_global_row_indices[0].size()!=f_VEC_ELEM.size()) {
      SETERRQ2(PETSC_COMM_SELF,1,"data inconsistency %d != %d",
	dof_global_row_indices[0].size(),f_VEC_ELEM.size());
    }
    ierr = VecSetValues(g,
      dof_global_row_indices[0].size(),&(dof_global_row_indices[0].data()[0]),
      &(f_VEC_ELEM.data())[0],ADD_VALUES); CHKERRQ(ierr);

  }
  PetscFunctionReturn(0);
}


}

