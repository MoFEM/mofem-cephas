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



#ifndef __PROJECTION10NODECOORDSONFIELD_HPP__
#define __PROJECTION10NODECOORDSONFIELD_HPP__

#include "FieldInterface.hpp"
#include "CoreDataStructures.hpp"

#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/matrix_proxy.hpp>
#include <boost/numeric/ublas/io.hpp>

using namespace boost::numeric;

namespace MoFEM {

struct Projection10NodeCoordsOnField: public FieldInterface::EntMethod { 

  FieldInterface& mField;
  string field_name;	
  int vErbose;

  Projection10NodeCoordsOnField(FieldInterface& _mField,string _field_name,int verb = 0): 
    mField(_mField),field_name(_field_name),vErbose(verb) {
  }

  PetscErrorCode ierr;
  ErrorCode rval;

  PetscErrorCode preProcess() {
    PetscFunctionBegin;
    PetscFunctionReturn(0);
  }

  ublas::vector<double> coords;
  ublas::vector<double,ublas::bounded_array<double,3> > ave_mid_coord;
  ublas::vector<double,ublas::bounded_array<double,3> > mid_node_coord;
  ublas::vector<double,ublas::bounded_array<double,3> > diff_node_coord;
  ublas::vector<double,ublas::bounded_array<double,3> > Dof;

  PetscErrorCode operator()() {
    PetscFunctionBegin;
    if(dof_ptr->get_name() != field_name) PetscFunctionReturn(0);
    if(dof_ptr->get_ent_type() == MBVERTEX) {
      EntityHandle node = dof_ptr->get_ent();
      coords.resize(3);
      rval = mField.get_moab().get_coords(&node,1,&*coords.data().begin());  CHKERR(rval);
      dof_ptr->get_FieldData() = coords[dof_ptr->get_dof_rank()];
      if(vErbose>0) {
	PetscPrintf(PETSC_COMM_WORLD,"val = %6.7e\n",dof_ptr->get_FieldData());
      }
      PetscFunctionReturn(0);
    }
    if(dof_ptr->get_ent_type() != MBEDGE) {
      PetscFunctionReturn(0);
    }
    if(dof_ptr->get_EntDofIdx() != dof_ptr->get_dof_rank()) {
      PetscFunctionReturn(0);
    }
    EntityHandle edge = dof_ptr->get_ent();
    if(mField.get_moab().type_from_handle(edge)!=MBEDGE) {
      SETERRQ(PETSC_COMM_SELF,1,"this method works only elements which are type of MBEDGE"); 
    }
    //coords
    int num_nodes;
    const EntityHandle* conn; 
    rval = mField.get_moab().get_connectivity(edge,conn,num_nodes,false); CHKERR_PETSC(rval);
    if( (num_nodes != 2) && (num_nodes != 3) ) {
      SETERRQ(PETSC_COMM_SELF,1,"this method works only 4 node and 10 node tets"); 
    }
    coords.resize(num_nodes*3);
    rval = mField.get_moab().get_coords(conn,num_nodes,&*coords.data().begin());  CHKERR(rval);
    ave_mid_coord.resize(3);
    mid_node_coord.resize(3);
    for(int dd = 0;dd<3;dd++) {
      ave_mid_coord[dd] = (coords[0*3+dd]+coords[1*3+dd])*0.5;
      if(num_nodes == 3) {
	mid_node_coord[dd] = coords[2*3+dd];
      } else {
	mid_node_coord[dd] = ave_mid_coord[dd];
      }
    }
    diff_node_coord.resize(3);
    ublas::noalias(diff_node_coord) = mid_node_coord-ave_mid_coord;
    double edge_shape_function_val = 0.25;
    // Dof*edge_shape_function_val + ave_mid_coord = mid_node_coord
    // Dof = (mid_node_coord-ave_mid_coord)/edge_shape_function_val
    Dof.resize(3);
    ublas::noalias(Dof) = diff_node_coord/edge_shape_function_val;
    /*if(dof_ptr->get_dof_order() != 2) {
      SETERRQ(PETSC_COMM_SELF,1,"this method works only fileds which are order 2"); 
    }*/
    if(dof_ptr->get_max_rank() != 3) {
      SETERRQ(PETSC_COMM_SELF,1,"this method works only fields which are rank 3"); 
    }
    dof_ptr->get_FieldData() = Dof[dof_ptr->get_dof_rank()];
    PetscFunctionReturn(0);
  }

  PetscErrorCode postProcess() {
    PetscFunctionBegin;
    PetscFunctionReturn(0);
  }

};

struct ProjectionFieldOn10NodeTet: public Projection10NodeCoordsOnField {

  bool set_nodes;
  bool on_coords;
  string on_tag;
  ProjectionFieldOn10NodeTet(FieldInterface& _mField,string _field_name,
    bool _set_nodes,bool _on_coords,string _on_tag = "NoNE"): 
    Projection10NodeCoordsOnField(_mField,_field_name),
    set_nodes(_set_nodes),on_coords(_on_coords),on_tag(_on_tag) {}
  
  Tag th;
  MoFEMField_multiIndex::index<FieldName_mi_tag>::type::iterator field_it;
  ublas::vector<double> L;
  PetscErrorCode preProcess() {
    PetscFunctionBegin;
    if(!on_coords) {
      if(on_tag == "NoNE") {
	SETERRQ(PETSC_COMM_SELF,1,"tag name not specified");
      }
      field_it = moabfields->get<FieldName_mi_tag>().find(field_name);
      if(field_it == moabfields->get<FieldName_mi_tag>().end()) {
	SETERRQ1(PETSC_COMM_SELF,1,"field not found %s",field_name.c_str());
      }
      int field_rank = field_it->get_max_rank();
      ublas::vector<double> def_VAL = ublas::zero_vector<double>(field_rank);
      rval = mField.get_moab().tag_get_handle(
	on_tag.c_str(),field_rank,MB_TYPE_DOUBLE,
	th,MB_TAG_CREAT|MB_TAG_SPARSE,&*def_VAL.data().begin()); CHKERR_THROW(rval);
    }
    L.resize(max_ApproximationOrder+1);
    ierr = Lagrange_basis(max_ApproximationOrder,0.,NULL,&*L.data().begin(),NULL,3); CHKERRQ(ierr);
    PetscFunctionReturn(0);
  }

  PetscErrorCode operator()() {
    PetscFunctionBegin;
    if(dof_ptr->get_name() != field_name) PetscFunctionReturn(0);
    if(set_nodes) {
      if(dof_ptr->get_ent_type() == MBVERTEX) {
	EntityHandle node = dof_ptr->get_ent();
	if(on_coords) {
	  coords.resize(3);
	  rval = mField.get_moab().get_coords(&node,1,&*coords.data().begin());  CHKERR(rval);
	  coords[dof_ptr->get_dof_rank()] = dof_ptr->get_FieldData();
	  rval = mField.get_moab().set_coords(&node,1,&*coords.data().begin());  CHKERR(rval);
	} else {
	  int field_rank = field_it->get_max_rank();
	  if(field_rank != dof_ptr->get_max_rank()) {
	    SETERRQ(PETSC_COMM_SELF,1,"data inconsistency");
	  }
	  double *tag_value;
	  int tag_size;
	  rval = mField.get_moab().tag_get_by_ptr(th,&node,1,(const void **)&tag_value,&tag_size); CHKERR_PETSC(rval);
	  if(tag_size != dof_ptr->get_max_rank()) {
	    SETERRQ(PETSC_COMM_SELF,1,"data inconsistency");
	  }
	  tag_value[dof_ptr->get_dof_rank()] = dof_ptr->get_FieldData();
	}
      }
      PetscFunctionReturn(0);
    }
    if(dof_ptr->get_ent_type() != MBEDGE) {
      PetscFunctionReturn(0);
    }
    EntityHandle edge = dof_ptr->get_ent();
    if(mField.get_moab().type_from_handle(edge)!=MBEDGE) {
      SETERRQ(PETSC_COMM_SELF,1,"this method works only elements which are type of MBEDGE"); 
    }

    int num_nodes;
    const EntityHandle* conn; 
    rval = mField.get_moab().get_connectivity(edge,conn,num_nodes,false); CHKERR_PETSC(rval);
    if( (num_nodes != 2) && (num_nodes != 3) ) {
      SETERRQ(PETSC_COMM_SELF,1,"this method works only 4 node and 10 node tets"); 
    }
    if(num_nodes == 2) {
      PetscFunctionReturn(0);
    }

    double approx_val = 0.25*L[dof_ptr->get_dof_order()-2]*dof_ptr->get_FieldData();;
    if(on_coords) {
      coords.resize(num_nodes*3);
      rval = mField.get_moab().get_coords(conn,num_nodes,&*coords.data().begin());  CHKERR(rval);
      if(dof_ptr->get_EntDofIdx() == dof_ptr->get_dof_rank()) { 
	//add only one when higher order terms present
	double ave_mid = (coords[3*0+dof_ptr->get_dof_rank()] + coords[3*1+dof_ptr->get_dof_rank()])*0.5;
	coords[2*3+dof_ptr->get_dof_rank()] = ave_mid;
      }
      coords[2*3+dof_ptr->get_dof_rank()] += approx_val;
      rval = mField.get_moab().set_coords(&conn[2],1,&coords[3*2]);  CHKERR(rval);
    } else {
      int tag_size;
      double *tag_value[num_nodes];
      rval = mField.get_moab().tag_get_by_ptr(
	th,conn,num_nodes,(const void **)tag_value,&tag_size); CHKERR_PETSC(rval);
      if(tag_size != dof_ptr->get_max_rank()) {
	SETERRQ(PETSC_COMM_SELF,1,"data inconsistency");
      }
      if(dof_ptr->get_EntDofIdx() == dof_ptr->get_dof_rank()) { 
	//add only one when higher order terms present
	double ave_mid = (tag_value[0][dof_ptr->get_dof_rank()] + tag_value[1][dof_ptr->get_dof_rank()])*0.5;
	tag_value[2][dof_ptr->get_dof_rank()] = ave_mid;
      }
      tag_value[2][dof_ptr->get_dof_rank()] += approx_val;
    }
    PetscFunctionReturn(0);
  }


};


}

#endif // __PROJECTION10NODECOORDSONFIELD_HPP__

