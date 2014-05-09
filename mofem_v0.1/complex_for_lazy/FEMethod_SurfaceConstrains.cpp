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
    ig_VEC_ELEM.resize(3);
    CT_MAT_ELEM.resize(9,3);
    dC_MAT_ELEM.resize(3,9);
    dCT_MAT_ELEM.resize(9,9);
    //crack front nodes
    Range crack_front_edges;
    mField.get_Cubit_msId_entities_by_dimension(201,SideSet,1,crack_front_edges,true);
    mField.get_moab().get_connectivity(crack_front_edges,crackFrontEdgesNodes,true);
    //projected nodes
    mField.get_moab().tag_get_handle(
      "PROJECTION_CRACK_SURFACE",3,MB_TYPE_DOUBLE,th_projection,MB_TAG_CREAT|MB_TAG_SPARSE,def_VAL);
  }
 
C_SURFACE_FEMethod::C_SURFACE_FEMethod(FieldInterface& _mField,BaseDirihletBC *_dirihlet_bc_method_ptr,Mat _C,string _lambda_field_name,int _verbose): 
    FEMethod(),mField(_mField),moab(_mField.get_moab()),dirihlet_bc_method_ptr(_dirihlet_bc_method_ptr),C(_C),lambda_field_name(_lambda_field_name) {
    run_in_constructor();
  }
C_SURFACE_FEMethod::C_SURFACE_FEMethod(FieldInterface& _mField,BaseDirihletBC *_dirihlet_bc_method_ptr,Mat _C,int _verbose): 
    FEMethod(),mField(_mField),moab(_mField.get_moab()),dirihlet_bc_method_ptr(_dirihlet_bc_method_ptr),C(_C),lambda_field_name("LAMBDA_SURFACE") {
    run_in_constructor();
  }

PetscErrorCode C_SURFACE_FEMethod::cOnstrain(double *dofs_iX,double *C,double *iC,double *g,double *ig) {
  PetscFunctionBegin;
  //set complex material position vector
  __CLPK_doublecomplex x_dofs_X0[9];
  __CLPK_doublecomplex x_dofs_X[9];
  for(int nn = 0;nn<3;nn++) {
    for(int dd = 0;dd<3;dd++) {
      if(lambda_global_row_indices[nn] == -1) {
	x_dofs_X0[nn*3+dd].r = ent_dofs_data[nn*3+dd];
      } else {
	x_dofs_X0[nn*3+dd].r = coords[nn*3+dd];
      }
      x_dofs_X[nn*3+dd].r = ent_dofs_data[nn*3+dd];
      if(dofs_iX != NULL) {
	if(lambda_global_row_indices[nn] == -1) {
	  x_dofs_X0[nn*3+dd].i = dofs_iX[nn*3+dd];
	} else {
	  x_dofs_X0[nn*3+dd].i = 0;
	}
	x_dofs_X[nn*3+dd].i = dofs_iX[nn*3+dd];
      } else {
	x_dofs_X0[nn*3+dd].i = 0;
	x_dofs_X[nn*3+dd].i = 0;
      }
  }}
  //calulate normal
  __CLPK_doublecomplex x_normal[3];
  ierr = ShapeFaceNormalMBTRI_complex(&diffNTRI[0],x_dofs_X,x_normal); CHKERRQ(ierr);
  //set direction if crack or interface surface
  Range adj_side_elems;
  BitRefLevel bit = problem_ptr->get_BitRefLevel();
  ierr = mField.get_adjacencies(bit,&face,1,3,adj_side_elems,Interface::INTERSECT,0); CHKERRQ(ierr);
  adj_side_elems = adj_side_elems.subset_by_type(MBTET);
  if(adj_side_elems.size()==0) {
    Range adj_tets_on_surface;
    BitRefLevel bit_tet_on_surface;
    bit_tet_on_surface.set(BITREFLEVEL_SIZE-1);
    ierr = mField.get_adjacencies(bit_tet_on_surface,&face,1,3,adj_tets_on_surface,Interface::INTERSECT,0); CHKERRQ(ierr);
    adj_side_elems.insert(*adj_tets_on_surface.begin());
  }
  if(adj_side_elems.size()!=1) {
    adj_side_elems.clear();
    ierr = mField.get_adjacencies(bit,&face,1,3,adj_side_elems,Interface::INTERSECT,5); CHKERRQ(ierr);
    Range::iterator it = adj_side_elems.begin();
    for(;it!=adj_side_elems.end();it++) {	
      Range nodes;
      rval = moab.get_connectivity(&*it,1,nodes,true); CHKERR_PETSC(rval);
      PetscPrintf(PETSC_COMM_WORLD,"%lu %lu %lu %lu\n",nodes[0],nodes[1],nodes[2],nodes[3]);
    }
    ParallelComm* pcomm = ParallelComm::get_pcomm(&mField.get_moab(),MYPCOMM_INDEX);
    if(pcomm->rank()==0) {
      EntityHandle out_meshset;
      rval = mField.get_moab().create_meshset(MESHSET_SET,out_meshset); CHKERR_PETSC(rval);
      rval = mField.get_moab().add_entities(out_meshset,adj_side_elems); CHKERR_PETSC(rval);
      rval = mField.get_moab().add_entities(out_meshset,&face,1); CHKERR_PETSC(rval);
      rval = mField.get_moab().write_file("debug_error.vtk","VTK","",&out_meshset,1); CHKERR_PETSC(rval);
    }
    SETERRQ1(PETSC_COMM_SELF,1,"expect 1 tet but is %u",adj_side_elems.size());
  }
  EntityHandle side_elem = *adj_side_elems.begin();
  if(side_elem!=0) {
    int side_number,sense,offset;
    rval = moab.side_number(side_elem,face,side_number,sense,offset); CHKERR_PETSC(rval);
    if(sense == -1) {
      __CLPK_doublecomplex xdot = { -1,0 };
      cblas_zscal(3,&xdot,x_normal,1);
    }
  }
  if(use_projection_from_crack_front) {
    Tag th_interface_side;
    rval = moab.tag_get_handle("INTERFACE_SIDE",th_interface_side); CHKERR_PETSC(rval);
    int side;
    rval = moab.tag_get_data(th_interface_side,&face,1,&side); CHKERR_PETSC(rval);
    if(side == 1) {
      __CLPK_doublecomplex xdot = { -1,0 };
      cblas_zscal(3,&xdot,x_normal,1);
    }
  }
  /*//save tag
  double normal0[3];
  ierr = ShapeFaceNormalMBTRI(&diffNTRI[0],&*coords.begin(),normal0); CHKERRQ(ierr);
  double area0 = cblas_dnrm2(3,normal0,1);
  cblas_dscal(3,1./area0,normal0,1);
  double def_NORMAL[9] = { 0,0,0 };
  Tag th_normal0;
  rval = moab.tag_get_handle("NORMAL0",3,MB_TYPE_DOUBLE,th_normal0,MB_TAG_CREAT|MB_TAG_SPARSE,def_NORMAL); CHKERR_PETSC(rval);
  if(side_elem!=0) {
    int side_number,sense,offset;
    rval = moab.side_number(side_elem,face,side_number,sense,offset); CHKERR_PETSC(rval);
    if(sense == -1) {
      cblas_dscal(3,-1,normal0,1);
    }
  }
  rval = moab.tag_set_data(th_normal0,&face,1,normal0); CHKERR_PETSC(rval);*/
  //calulare complex normal length
  /*double __complex__ xarea = csqrt(
      cpow((x_normal[0].r+I*x_normal[0].i),2)+
      cpow((x_normal[1].r+I*x_normal[1].i),2)+
      cpow((x_normal[2].r+I*x_normal[2].i),2));*/
  //
  if( C!=NULL) bzero( C,3*9*sizeof(double));
  if(iC!=NULL) bzero(iC,3*9*sizeof(double));
  if( g!=NULL) bzero( g,3*sizeof(double));
  if(ig!=NULL) bzero(ig,3*sizeof(double));
  for(unsigned int gg = 0;gg<g_NTRI.size()/3;gg++) {
    for(int nn = 0;nn<3;nn++) {
      for(int mm = 0;mm<3;mm++) {
	for(int dd = 0;dd<3;dd++) {
	  double __complex__ c = (x_normal [dd].r+I*x_normal [dd].i)*g_NTRI[3*gg+mm];
	  double __complex__ xg  = c*(x_dofs_X [3*mm+dd].r+I*x_dofs_X [3*mm+dd].i);
	  double __complex__ xg0 = c*(x_dofs_X0[3*mm+dd].r+I*x_dofs_X0[3*mm+dd].i);
	  if( C!=NULL) {
	     C[9*nn + 3*mm+dd] += G_TRI_W[gg]*g_NTRI[3*gg+nn]*creal(c);
	  }
	  if(iC!=NULL) {
	    iC[9*nn + 3*mm+dd] += G_TRI_W[gg]*g_NTRI[3*gg+nn]*cimag(c);
	  }
	  if( g!=NULL) {
	     g[nn] +=  G_TRI_W[gg]*g_NTRI[3*gg+nn]*creal(xg-xg0);
	  }
	  if(ig!=NULL) {
	    ig[nn] +=  G_TRI_W[gg]*g_NTRI[3*gg+nn]*cimag(xg-xg0);
	  }
	}
      }
    }
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
    ierr = cOnstrain(NULL,&*C_MAT_ELEM.data().begin(),NULL,NULL,NULL); CHKERRQ(ierr);
    if(transpose) {
      ublas::noalias(CT_MAT_ELEM) = trans(C_MAT_ELEM);
    } 
    if(nonlinear) {
      ublas::noalias(dC_MAT_ELEM) = ublas::zero_matrix<double>(3,9);
      ublas::noalias(dCT_MAT_ELEM) = ublas::zero_matrix<double>(9,9);
      for(int dd = 0;dd<9;dd++) {
	ublas::noalias(ent_idofs_data) = ublas::zero_vector<double>(9);
	ent_idofs_data[dd] = r*eps;
	ierr = cOnstrain(&*ent_idofs_data.data().begin(),NULL,&*iC_MAT_ELEM.data().begin(),NULL,&*ig_VEC_ELEM.data().begin()); CHKERRQ(ierr);
	for(int nnn = 0;nnn<3;nnn++) {
	  dC_MAT_ELEM(nnn,dd) = ig_VEC_ELEM[nnn]/(r*eps);
	}
      	if_VEC_ELEM = prod(trans(iC_MAT_ELEM)/(r*eps),ent_lambda_data);
	for(int nnn = 0;nnn<3;nnn++) {
	  for(int ddd = 0;ddd<3;ddd++) {
	    dCT_MAT_ELEM(3*nnn+ddd,dd) = if_VEC_ELEM[3*nnn+ddd];
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
    if(nonlinear) {
      ierr = MatSetValues(C,
	lambda_global_row_indices.size(),&(lambda_global_row_indices.data()[0]),
	dof_global_col_indices[0].size(),&(dof_global_col_indices[0].data()[0]),
	&(dC_MAT_ELEM.data())[0],ADD_VALUES); CHKERRQ(ierr);
    } else {
      ierr = MatSetValues(C,
	lambda_global_row_indices.size(),&(lambda_global_row_indices.data()[0]),
	dof_global_col_indices[0].size(),&(dof_global_col_indices[0].data()[0]),
	&(C_MAT_ELEM.data())[0],ADD_VALUES); CHKERRQ(ierr);
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
	  if(lambda_global_row_indices[nn] == -1) {
	    dof_global_row_indices[0][nn*3+0] = -1;
	    dof_global_row_indices[0][nn*3+1] = -1;
	    dof_global_row_indices[0][nn*3+2] = -1;
	  } else {
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
	}
	DirihletBC.clear();
	ierr = dirihlet_bc_method_ptr->SetDirihletBC_to_ElementIndicies(
	  this,dof_global_row_indices,dof_global_col_indices,DirihletBC); CHKERRQ(ierr);
      }
      rval = moab.get_coords(conn_face,num_nodes,&*coords.data().begin()); CHKERR_PETSC(rval);
      if(use_projection_from_crack_front) {
	for(int nn = 0;nn<3;nn++) {
	  if(crackFrontEdgesNodes.find(conn_face[nn])!=crackFrontEdgesNodes.end()) {
	  double projection[3];
	    rval = mField.get_moab().tag_get_data(th_projection,&conn_face[nn],1,projection); CHKERR_PETSC(rval);
	    for(int dd = 0;dd<3;dd++) {
	      coords[nn*3+dd] = projection[dd];
	    }
	  }
	}
      }
      if(!nonlinear) {
	ublas::noalias(ent_dofs_data) = coords;
      }
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

g_SURFACE_FEMethod::g_SURFACE_FEMethod(FieldInterface& _mField,BaseDirihletBC *_dirihlet_bc_method_ptr,Vec _g,string _lambda_field_name,int _verbose): 
    C_SURFACE_FEMethod(_mField,_dirihlet_bc_method_ptr,PETSC_NULL,_lambda_field_name,_verbose),g(_g) {
    g_VEC_ELEM.resize(3);
  }
g_SURFACE_FEMethod::g_SURFACE_FEMethod(FieldInterface& _mField,BaseDirihletBC *_dirihlet_bc_method_ptr,Vec _g,int _verbose): 
    C_SURFACE_FEMethod(_mField,_dirihlet_bc_method_ptr,PETSC_NULL,_verbose),g(_g) {
    g_VEC_ELEM.resize(3);
  }

PetscErrorCode g_SURFACE_FEMethod::iNtegrate(bool transpose,bool nonlinear) {
  PetscFunctionBegin;
  try {
    if(transpose) {
      ierr = cOnstrain(NULL,&*C_MAT_ELEM.data().begin(),NULL,g_VEC_ELEM.data().begin(),NULL); CHKERRQ(ierr);
      CT_MAT_ELEM = trans(C_MAT_ELEM);
      f_VEC_ELEM = prod(CT_MAT_ELEM,ent_lambda_data);
    } else {
      ierr = cOnstrain(NULL,NULL,NULL,g_VEC_ELEM.data().begin(),NULL); CHKERRQ(ierr);
    }
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

