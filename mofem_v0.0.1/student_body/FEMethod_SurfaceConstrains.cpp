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

namespace MoFEM {

const int debug_constrains = 1;

PetscErrorCode C_SURFACE_FEMethod::preProcess() {
    PetscFunctionBegin;
    PetscFunctionReturn(0);
  }

PetscErrorCode C_SURFACE_FEMethod::SaveConstrainOnTags() {
    PetscFunctionBegin;
    PetscFunctionReturn(0);
  }

void C_SURFACE_FEMethod::run_in_constructor() {
    diffNTRI.resize(6);
    ShapeDiffMBTRI(&diffNTRI[0]);
    g_NTRI3.resize(4*3);
    ShapeMBTRI(&g_NTRI3[0],G_TRI_X4,G_TRI_Y4,4);
    G_TRI_W = G_TRI_W4;
    double def_VAL[3*9];
    fill(&def_VAL[0],&def_VAL[3*9],0);
    rval = moab.tag_get_handle("MATERIAL_NORMAL",3,MB_TYPE_DOUBLE,th_material_normal,MB_TAG_CREAT|MB_TAG_SPARSE,&def_VAL); CHKERR_THROW(rval);
    ent_global_col_indices.resize(9);
    ent_global_row_indices.resize(3);
    ent_lambda_data.resize(3);
    ent_dofs_data.resize(9);
    coords.resize(9);
    ent_normal_map0.resize(3);
    ent_normal_map.resize(3);
    C_MAT_ELEM.resize(3,9);
    CT_MAT_ELEM.resize(9,3);
  }
 
C_SURFACE_FEMethod::C_SURFACE_FEMethod(Interface& _moab,Mat _C,string _lambda_field_name,int _verbose): 
    FEMethod(),moab(_moab),C(_C),lambda_field_name(_lambda_field_name) {
    run_in_constructor();
  }
C_SURFACE_FEMethod::C_SURFACE_FEMethod(Interface& _moab,Mat _C,int _verbose): 
    FEMethod(),moab(_moab),C(_C),lambda_field_name("LAMBDA_SURFACE") {
    run_in_constructor();
  }

PetscErrorCode C_SURFACE_FEMethod::Integrate(bool transpose) {
  PetscFunctionBegin;
  try {
    ublas::noalias(C_MAT_ELEM) = ublas::zero_matrix<double>(3,9);
    double area0 = norm_2(ent_normal_map0);
    double area = norm_2(ent_normal_map);
    for(unsigned int gg = 0;gg<g_NTRI3.size()/3;gg++) {
	for(int nn = 0;nn<3;nn++) {
	  for(int dd = 0;dd<3;dd++) {
	    for(int nnn = 0;nnn<3;nnn++) {
	      C_MAT_ELEM(nn,3*nnn+dd) += 
		G_TRI_W[gg]*
		g_NTRI3[3*gg+nn]*g_NTRI3[3*gg+nnn]*
		(area0/area)*ent_normal_map[dd];
	    }
	  }
	}
    }
    if(transpose) {
      ublas::noalias(CT_MAT_ELEM) = trans(C_MAT_ELEM);
    }
    ierr = SaveConstrainOnTags(); CHKERRQ(ierr);
  } catch (const std::exception& ex) {
    ostringstream ss;
    ss << "thorw in method: " << ex.what() << " at line " << __LINE__ << " in file " << __FILE__ << endl;
    SETERRQ(PETSC_COMM_SELF,1,ss.str().c_str());
  }
  PetscFunctionReturn(0);
}

PetscErrorCode C_SURFACE_FEMethod::Assemble(bool transpose) {
    PetscFunctionBegin;
    ierr = MatSetValues(C,
      ent_global_row_indices.size(),&(ent_global_row_indices.data()[0]),
      ent_global_col_indices.size(),&(ent_global_col_indices.data()[0]),
      &(C_MAT_ELEM.data())[0],ADD_VALUES); CHKERRQ(ierr);
    if(transpose) {
      ierr = MatSetValues(C,
	ent_global_col_indices.size(),&(ent_global_col_indices.data()[0]),
	ent_global_row_indices.size(),&(ent_global_row_indices.data()[0]),
	&CT_MAT_ELEM.data()[0],ADD_VALUES); CHKERRQ(ierr);
    }
    PetscFunctionReturn(0);
}

PetscErrorCode C_SURFACE_FEMethod::operator()(bool transpose) {
    PetscFunctionBegin;
    try {
      EntityHandle face = fe_ptr->get_ent();
      const EntityHandle* conn_face; 
      int num_nodes; 
      rval = moab.get_connectivity(face,conn_face,num_nodes,true); CHKERR_PETSC(rval);
      for(int nn = 0;nn<num_nodes;nn++) {
  	FENumeredDofMoFEMEntity_multiIndex::index<Composite_mi_tag3>::type::iterator dit,hi_dit;
  	dit = row_multiIndex->get<Composite_mi_tag3>().lower_bound(boost::make_tuple(lambda_field_name,conn_face[nn]));
  	hi_dit = row_multiIndex->get<Composite_mi_tag3>().upper_bound(boost::make_tuple(lambda_field_name,conn_face[nn]));
  	if(distance(dit,hi_dit)==0) {
  	  ent_global_row_indices[nn] = -1;
  	} else {
  	  if(distance(dit,hi_dit)!=1) SETERRQ1(PETSC_COMM_SELF,1,"data inconsistency, number of dof on node for < %s > should be 1",lambda_field_name.c_str());
  	  int global_idx = dit->get_petsc_gloabl_dof_idx();
  	  ent_global_row_indices[nn] = global_idx;
  	  ent_lambda_data[nn] = dit->get_FieldData();
  	  int local_idx = dit->get_petsc_local_dof_idx();
  	  if(local_idx<0) SETERRQ(PETSC_COMM_SELF,1,"data inconsistency, negative index of local dofs on element");
  	}
  	dit = col_multiIndex->get<Composite_mi_tag3>().lower_bound(boost::make_tuple("MESH_NODE_POSITIONS",conn_face[nn]));
  	hi_dit = col_multiIndex->get<Composite_mi_tag3>().upper_bound(boost::make_tuple("MESH_NODE_POSITIONS",conn_face[nn]));
  	if(distance(dit,hi_dit)!=3) SETERRQ(PETSC_COMM_SELF,1,"data inconsistency, number of dof on node for MESH_NODE_POSITIONS should be 3");
  	  for(;dit!=hi_dit;dit++) {
  	    int global_idx = dit->get_petsc_gloabl_dof_idx();
  	    assert(nn*3+dit->get_dof_rank()<9);
  	    if(ent_global_row_indices[nn] == -1) {
  	      ent_global_col_indices[nn*3+dit->get_dof_rank()] = -1;
  	    } else {
  	      ent_global_col_indices[nn*3+dit->get_dof_rank()] = global_idx;
  	    }
  	    int local_idx = dit->get_petsc_local_dof_idx();
  	    if(local_idx<0) SETERRQ(PETSC_COMM_SELF,1,"data inconsistency, negative index of local dofs on element");
  	    ent_dofs_data[nn*3+dit->get_dof_rank()] = dit->get_FieldData();
  	}
      }
      rval = moab.get_coords(conn_face,num_nodes,&*coords.data().begin()); CHKERR_PETSC(rval);
      ierr = ShapeFaceNormalMBTRI(&diffNTRI[0],&coords.data()[0],&ent_normal_map0.data()[0]); CHKERRQ(ierr);
      ierr = ShapeFaceNormalMBTRI(&diffNTRI[0],&coords.data()[0],&ent_normal_map.data()[0]); CHKERRQ(ierr);
      ublas::vector<double> scaled_normal = ent_normal_map0/norm_2(ent_normal_map0);
      rval = moab.tag_set_data(th_material_normal,&face,1,&scaled_normal.data()[0]); CHKERR_PETSC(rval);
      ierr = this->Integrate(transpose); CHKERRQ(ierr);
      ierr = this->Assemble(transpose); CHKERRQ(ierr);
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

g_SURFACE_FEMethod::g_SURFACE_FEMethod(Interface& _moab,Vec _g,string _lambda_field_name,int _verbose): 
    C_SURFACE_FEMethod(_moab,PETSC_NULL,_lambda_field_name,_verbose),g(_g) {
    g_VEC_ELEM.resize(3);
  }
g_SURFACE_FEMethod::g_SURFACE_FEMethod(Interface& _moab,Vec _g,int _verbose): 
    C_SURFACE_FEMethod(_moab,PETSC_NULL,_verbose),g(_g) {
    g_VEC_ELEM.resize(3);
  }

PetscErrorCode g_SURFACE_FEMethod::Integrate(bool transpose) {
    PetscFunctionBegin;
    try {
      C_SURFACE_FEMethod::Integrate(transpose);
      /*ublas::noalias(g_VEC_ELEM) = ublas::zero_vector<double>(3);
      double area0 = norm_2(ent_normal_map0);
      double area = norm_2(ent_normal_map);
      for(unsigned int gg = 0;gg<g_NTRI3.size()/3;gg++) {
  	for(int nn = 0;nn<3;nn++) {
  	  for(int dd = 0;dd<3;dd++) {
  	    double X0_dd = cblas_ddot(3,&g_NTRI3[3*gg],1,&coords.data()[dd],3);
  	    double X_dd = cblas_ddot(3,&g_NTRI3[3*gg],1,&ent_dofs_data.data()[dd],3);
  	    g_VEC_ELEM[nn] += G_TRI_W[gg]*
		g_NTRI3[3*gg+nn]*
		  (ent_normal_map[dd]*X_dd*(area0/area) - ent_normal_map0[dd]*X0_dd);
  	  }
  	}
      }*/
      g_VEC_ELEM = prod(C_MAT_ELEM,ent_dofs_data-coords);
      f_VEC_ELEM = prod(CT_MAT_ELEM,ent_lambda_data);
    } catch (const std::exception& ex) {
      ostringstream ss;
      ss << "thorw in method: " << ex.what() << " at line " << __LINE__ << " in file " << __FILE__ << endl;
      SETERRQ(PETSC_COMM_SELF,1,ss.str().c_str());
    }
    PetscFunctionReturn(0);
}
PetscErrorCode g_SURFACE_FEMethod::Assemble(bool transpose) {
  PetscFunctionBegin;
  if(ent_global_row_indices.size()!=g_VEC_ELEM.size()) {
    SETERRQ2(PETSC_COMM_SELF,1,"data inconsistency %d != %d",
	ent_global_row_indices.size(),g_VEC_ELEM.size());
  }
  ierr = VecSetOption(g,VEC_IGNORE_NEGATIVE_INDICES,PETSC_TRUE); CHKERRQ(ierr);
  ierr = VecSetValues(g,
        ent_global_row_indices.size(),&(ent_global_row_indices.data()[0]),
        &(g_VEC_ELEM.data())[0],ADD_VALUES); CHKERRQ(ierr);
  if(transpose) {
    if(ent_global_col_indices.size()!=f_VEC_ELEM.size()) {
      SETERRQ2(PETSC_COMM_SELF,1,"data inconsistency %d != %d",
	ent_global_col_indices.size(),f_VEC_ELEM.size());
    }
    ierr = VecSetValues(g,
	ent_global_col_indices.size(),&(ent_global_col_indices.data()[0]),
	&(f_VEC_ELEM.data())[0],ADD_VALUES); CHKERRQ(ierr);

  }
  PetscFunctionReturn(0);
}

C_CORNER_FEMethod::C_CORNER_FEMethod(Interface& _moab,Mat _C,int _verbose): 
    FEMethod(),moab(_moab),C(_C) {
    C_MAT_ELEM.resize(3,3);
    ent_lambda_data.resize(3);
    ent_dofs_data.resize(3);
    ent_global_col_indices.resize(3);
    ent_global_row_indices.resize(3);
    coords.resize(3);
}

PetscErrorCode C_CORNER_FEMethod::preProcess() {
    PetscFunctionBegin;
    PetscFunctionReturn(0);
  }

PetscErrorCode C_CORNER_FEMethod::Integrate(bool transpose) {
    PetscFunctionBegin;
    try {
    ublas::noalias(C_MAT_ELEM) = ublas::zero_matrix<double>(3,3);
    for(int nn = 0;nn<3;nn++) {
	for(int dd = 0;dd<3;dd++) {
	  if(nn!=dd) C_MAT_ELEM(nn,dd) = 0;
	  else C_MAT_ELEM(nn,dd) = 1;
	}
    }
    ierr = MatSetValues(C,
      ent_global_row_indices.size(),&(ent_global_row_indices[0]),
      ent_global_col_indices.size(),&(ent_global_col_indices[0]),
      &(C_MAT_ELEM.data())[0],ADD_VALUES); CHKERRQ(ierr);
    if(transpose) {
      ierr = MatSetValues(C,
	ent_global_col_indices.size(),&(ent_global_col_indices[0]),
	ent_global_row_indices.size(),&(ent_global_row_indices[0]),
	&(C_MAT_ELEM.data())[0],ADD_VALUES); CHKERRQ(ierr);
    }
    } catch (const std::exception& ex) {
      ostringstream ss;
      ss << "thorw in method: " << ex.what() << " at line " << __LINE__ << " in file " << __FILE__ << endl;
      SETERRQ(PETSC_COMM_SELF,1,ss.str().c_str());
    }
    PetscFunctionReturn(0);
  }

PetscErrorCode C_CORNER_FEMethod::operator()(bool transpose) {
    PetscFunctionBegin;
    try {
    EntityHandle node = fe_ptr->get_ent();
    FENumeredDofMoFEMEntity_multiIndex::index<Composite_mi_tag3>::type::iterator dit,hi_dit;
    dit = col_multiIndex->get<Composite_mi_tag3>().lower_bound(boost::make_tuple("MESH_NODE_POSITIONS",node));
    hi_dit = col_multiIndex->get<Composite_mi_tag3>().upper_bound(boost::make_tuple("MESH_NODE_POSITIONS",node));
    if(distance(dit,hi_dit)!=3) {
      SETERRQ1(PETSC_COMM_SELF,1,"data inconsistency, number of dof on node for MESH_NODE_POSITIONS should be 3, but is %u",distance(dit,hi_dit));
    }
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
	ent_lambda_data[dit->get_dof_rank()] = dit->get_FieldData();
	if(local_idx<0) SETERRQ(PETSC_COMM_SELF,1,"data inconsistency, negative index of local dofs on element");
    }
    rval = moab.get_coords(&node,1,&*coords.data().begin()); CHKERR_PETSC(rval);
    ierr = this->Integrate(transpose); CHKERRQ(ierr);
    } catch (const std::exception& ex) {
      ostringstream ss;
      ss << "thorw in method: " << ex.what() << " at line " << __LINE__ << " in file " << __FILE__ << endl;
      SETERRQ(PETSC_COMM_SELF,1,ss.str().c_str());
    }
    PetscFunctionReturn(0);
  }

PetscErrorCode C_CORNER_FEMethod::postProcess() {
    PetscFunctionBegin;
    PetscFunctionReturn(0);
  }

g_CORNER_FEMethod::g_CORNER_FEMethod(Interface& _moab,Vec _g,int _verbose): 
    C_CORNER_FEMethod(_moab,PETSC_NULL,_verbose),g(_g) {
    g_VEC_ELEM.resize(3);
  }

PetscErrorCode g_CORNER_FEMethod::Integrate(bool transpose) {
    PetscFunctionBegin;
    try {
      for(int nn = 0;nn<3;nn++) {
        g_VEC_ELEM[nn] = ent_dofs_data[nn] - coords[nn];
      }
      VecSetOption(g,VEC_IGNORE_NEGATIVE_INDICES,PETSC_TRUE); 
      ierr = VecSetValues(g,
        ent_global_row_indices.size(),&(ent_global_row_indices[0]),
        &(g_VEC_ELEM.data())[0],ADD_VALUES); CHKERRQ(ierr);
      if(transpose) {
	if(ent_lambda_data.size()!=ent_global_col_indices.size()) {
	  SETERRQ(PETSC_COMM_SELF,1,"data inconsistency");
	}
	ierr = VecSetValues(g,
	  ent_global_col_indices.size(),&(ent_global_col_indices[0]),
	  &(ent_lambda_data.data())[0],ADD_VALUES); CHKERRQ(ierr);
      }
    } catch (const std::exception& ex) {
      ostringstream ss;
      ss << "thorw in method: " << ex.what() << " at line " << __LINE__ << " in file " << __FILE__ << endl;
      SETERRQ(PETSC_COMM_SELF,1,ss.str().c_str());
    }
    PetscFunctionReturn(0);
  }

}

