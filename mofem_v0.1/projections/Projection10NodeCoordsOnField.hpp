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

struct Projection10NodeCoordsOnField: public FieldInterface::EntMethod { //FieldInterface::FEMethod {

  FieldInterface& mField;
  string field_name;

  Projection10NodeCoordsOnField(FieldInterface& _mField,string _field_name): 
    mField(_mField),field_name(_field_name) {
  }

  PetscErrorCode ierr;
  ErrorCode rval;

  Range test;
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
      PetscFunctionReturn(0);
    }
    if(dof_ptr->get_ent_type() != MBEDGE) {
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
    double edge_shape_function_val = 0.5*0.5;
    // Dof*edge_shape_function_val + ave_mid_coord = mid_node_coord
    // Dof = (mid_node_coord-ave_mid_coord)/edge_shape_function_val
    Dof.resize(3);
    ublas::noalias(Dof) = diff_node_coord/edge_shape_function_val;
    if(dof_ptr->get_dof_order() != 2) {
      SETERRQ(PETSC_COMM_SELF,1,"this method works only fileds which are order 2"); 
    }
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

}

#endif // __PROJECTION10NODECOORDSONFIELD_HPP__

