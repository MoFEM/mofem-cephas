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


struct SurfaceFEMethod:public moabField::FEMethod {
  ErrorCode rval;
  PetscErrorCode ierr;
  Interface& moab;

  Tag th_material_normal;
  Tag th_const_shape;
  vector<double> diffNTRI;
  vector<double> g_NTRI3;
  SurfaceFEMethod(Interface& _moab,int _verbose = 0): FEMethod(),moab(_moab) {
    diffNTRI.resize(6);
    ShapeDiffMBTRI(&diffNTRI[0]);
    g_NTRI3.resize(3*3);
    ShapeMBTRI(&g_NTRI3[0],G_TRI_X3,G_TRI_Y3,3);
    double def_VAL[3] = { 0,0,0 };
    rval = moab.tag_get_handle("MATERIAL_NORMAL",3,MB_TYPE_DOUBLE,th_material_normal,MB_TAG_CREAT|MB_TAG_SPARSE,&def_VAL); CHKERR_THROW(rval);
    rval = moab.tag_get_handle("CONST_SHAPE",3,MB_TYPE_DOUBLE,th_const_shape,MB_TAG_CREAT|MB_TAG_SPARSE,&def_VAL); CHKERR_THROW(rval);
  }

  PetscErrorCode preProcess() {
    PetscFunctionBegin;
    PetscFunctionReturn(0);
  }

  ublas::vector<double> C_CONST_SHAPE;
  Range skin_faces;
  map<EntityHandle,ublas::vector<DofIdx,ublas::bounded_array<DofIdx,9> > > ent_global_col_indices;
  map<EntityHandle,ublas::vector<double,ublas::bounded_array<double,9> > > ent_dofs_data;
  map<EntityHandle,ublas::vector<double,ublas::bounded_array<double,3> > > ent_normal_map;
  PetscErrorCode operator()() {
    PetscFunctionBegin;
    SideNumber_multiIndex &side_table = fe_ptr->get_side_number_table();
    SideNumber_multiIndex::nth_index<1>::type::iterator siit = side_table.get<1>().lower_bound(boost::make_tuple(MBTRI,0));
    SideNumber_multiIndex::nth_index<1>::type::iterator hi_siit = side_table.get<1>().upper_bound(boost::make_tuple(MBTRI,4));
    ent_normal_map.clear();
    ent_global_col_indices.clear();
    ent_dofs_data.clear();
    for(;siit!=hi_siit;siit++) {
      EntityHandle face = siit->ent;
      Range adj_tets;
      rval = moab.get_adjacencies(&face,1,3,false,adj_tets); CHKERR_PETSC(rval);
      if(adj_tets.size()>1) continue;
      ent_dofs_data[face].resize(9);
      ent_global_col_indices[face].resize(9);
      const EntityHandle* conn_face; 
      int num_nodes; 
      rval = moab.get_connectivity(face,conn_face,num_nodes,true); CHKERR_PETSC(rval);
      for(int nn = 0;nn<num_nodes;nn++) {
	FENumeredDofMoFEMEntity_multiIndex::index<Composite_mi_tag3>::type::iterator dit,hi_dit;
	dit = col_multiIndex->get<Composite_mi_tag3>().lower_bound(boost::make_tuple("MESH_NODE_POSITIONS",conn_face[nn]));
	hi_dit = col_multiIndex->get<Composite_mi_tag3>().upper_bound(boost::make_tuple("MESH_NODE_POSITIONS",conn_face[nn]));
	for(;dit!=hi_dit;dit++) {
	  ent_dofs_data[face][nn*3+dit->get_dof_rank()] = dit->get_FieldData();
	  ent_global_col_indices[face][nn*3+dit->get_dof_rank()] = dit->get_petsc_gloabl_dof_idx();
	}
      }
      ent_normal_map[face].resize(3);
      ierr = ShapeFaceNormalMBTRI(&diffNTRI[0],&ent_dofs_data[face].data()[0],&ent_normal_map[face].data()[0]); CHKERRQ(ierr);
      rval = moab.tag_set_data(th_material_normal,&face,1,&ent_normal_map[face].data()[0]); CHKERR_PETSC(rval);
      skin_faces.insert(face);
      C_CONST_SHAPE.resize(3*3);
      for(int nn = 0;nn<num_nodes;nn++) {
	for(int dd = 0;dd<3;dd++) {
	  C_CONST_SHAPE[nn*3+dd] = ent_normal_map[face][dd];
	}
      }
      C_CONST_SHAPE /= norm_2(ent_normal_map[face]);
    }
    PetscFunctionReturn(0);
  }
  PetscErrorCode postProcess() {
    PetscFunctionBegin;
    PetscFunctionReturn(0);
  }

};

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
