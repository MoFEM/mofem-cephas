/** \file Projection10NodeCoordsOnField.hpp

FIXME: Move code to cpp file.

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

using namespace boost::numeric;

namespace MoFEM {

/** \brief Projection of edge entities with one mid-node on hierarchical basis
*/
struct Projection10NodeCoordsOnField: public EntMethod {

  Interface& mField;
  std::string field_name;
  int vErbose;

  Projection10NodeCoordsOnField(MoFEM::Interface& m_field,std::string _field_name,int verb = 0):
    mField(m_field),field_name(_field_name),vErbose(verb) {
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
    if(dofPtr == NULL) {
      SETERRQ(PETSC_COMM_SELF,MOFEM_DATA_INCONSISTENCY,"data inconsistency");
    }
    if(dofPtr->getName() != field_name) PetscFunctionReturn(0);
    if(dofPtr->getEntType() == MBVERTEX) {
      EntityHandle node = dofPtr->getEnt();
      coords.resize(3);
      rval = mField.get_moab().get_coords(&node,1,&*coords.data().begin());  CHKERR_MOAB(rval);
      dofPtr->getFieldData() = coords[dofPtr->getDofCoeffIdx()];
      if(vErbose>0) {
        PetscPrintf(mField.get_comm(),"val = %6.7e\n",dofPtr->getFieldData());
      }
      PetscFunctionReturn(0);
    }
    if(dofPtr->getEntType() != MBEDGE) {
      PetscFunctionReturn(0);
    }
    if(dofPtr->getEntDofIdx() != dofPtr->getDofCoeffIdx()) {
      PetscFunctionReturn(0);
    }
    EntityHandle edge = dofPtr->getEnt();
    if(mField.get_moab().type_from_handle(edge)!=MBEDGE) {
      SETERRQ(PETSC_COMM_SELF,1,"this method works only elements which are type of MBEDGE");
    }
    //coords
    int num_nodes;
    const EntityHandle* conn;
    rval = mField.get_moab().get_connectivity(edge,conn,num_nodes,false); CHKERRQ_MOAB(rval);
    if( (num_nodes != 2) && (num_nodes != 3) ) {
      SETERRQ(PETSC_COMM_SELF,1,"this method works only 4 node and 10 node tets");
    }
    coords.resize(num_nodes*3);
    rval = mField.get_moab().get_coords(conn,num_nodes,&*coords.data().begin());  CHKERR_MOAB(rval);
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
    double edge_shape_function_val = 0.25;
    FieldApproximationBase base = dofPtr->getApproxBase();
    switch (base) {
      case AINSWORTH_COLE_BASE:
      break;
      case LOBATTO_BASE:
      edge_shape_function_val *= LOBATTO_PHI0(0);
      break;
      default:
      SETERRQ(PETSC_COMM_SELF,MOFEM_NOT_IMPLEMENTED,"not yet implemented");
    }

    diff_node_coord.resize(3);
    ublas::noalias(diff_node_coord) = mid_node_coord-ave_mid_coord;
    // Dof*edge_shape_function_val + ave_mid_coord = mid_node_coord
    // Dof = (mid_node_coord-ave_mid_coord)/edge_shape_function_val
    Dof.resize(3);
    ublas::noalias(Dof) = diff_node_coord/edge_shape_function_val;
    /*if(dofPtr->getDofOrder() != 2) {
      SETERRQ(PETSC_COMM_SELF,1,"this method works only fileds which are order 2");
    }*/
    if(dofPtr->getNbOfCoeffs() != 3) {
      SETERRQ(PETSC_COMM_SELF,1,"this method works only fields which are rank 3");
    }
    dofPtr->getFieldData() = Dof[dofPtr->getDofCoeffIdx()];
    PetscFunctionReturn(0);
  }

  PetscErrorCode postProcess() {
    PetscFunctionBegin;
    PetscFunctionReturn(0);
  }

};

struct ProjectionFieldOn10NodeTet: public Projection10NodeCoordsOnField {

  bool setNodes;
  bool onCoords;
  std::string onTag;

  const int maxApproximationOrder;

  ProjectionFieldOn10NodeTet(
    MoFEM::Interface& m_field,std::string _field_name,bool set_nodes,bool on_coords,std::string on_tag = "NoNE"
  ):
  Projection10NodeCoordsOnField(m_field,_field_name),
  setNodes(set_nodes),
  onCoords(on_coords),
  onTag(on_tag),
  maxApproximationOrder(20) {
  }

  Tag th;
  Field_multiIndex::index<FieldName_mi_tag>::type::iterator field_it;
  ublas::vector<double> L;
  ublas::vector<double> K;

  PetscErrorCode preProcess() {
    PetscFunctionBegin;
    if(!onCoords) {
      if(onTag == "NoNE") {
        SETERRQ(PETSC_COMM_SELF,1,"tag name not specified");
      }
      field_it = fieldsPtr->get<FieldName_mi_tag>().find(field_name);
      if(field_it == fieldsPtr->get<FieldName_mi_tag>().end()) {
        SETERRQ1(PETSC_COMM_SELF,1,"field not found %s",field_name.c_str());
      }
      int field_rank = (*field_it)->getNbOfCoeffs();
      ublas::vector<double> def_VAL = ublas::zero_vector<double>(field_rank);
      rval = mField.get_moab().tag_get_handle(
        onTag.c_str(),field_rank,MB_TYPE_DOUBLE,
        th,MB_TAG_CREAT|MB_TAG_SPARSE,&*def_VAL.data().begin()
      ); MOAB_THROW(rval);
    }

    L.resize(maxApproximationOrder+1);
    ierr = Legendre_polynomials(maxApproximationOrder,0.,NULL,&*L.data().begin(),NULL,3); CHKERRQ(ierr);
    K.resize(10);
    ierr = LobattoKernel_polynomials(9,0.,NULL,&*K.data().begin(),NULL,3); CHKERRQ(ierr);

    PetscFunctionReturn(0);
  }

  PetscErrorCode operator()() {
    PetscFunctionBegin;
    if(dofPtr->getName() != field_name) PetscFunctionReturn(0);
    if(setNodes) {
      if(dofPtr->getEntType() == MBVERTEX) {
        EntityHandle node = dofPtr->getEnt();
        if(onCoords) {
          coords.resize(3);
          rval = mField.get_moab().get_coords(&node,1,&*coords.data().begin());  CHKERR_MOAB(rval);
          coords[dofPtr->getDofCoeffIdx()] = dofPtr->getFieldData();
          rval = mField.get_moab().set_coords(&node,1,&*coords.data().begin());  CHKERR_MOAB(rval);
        } else {
          int field_rank = (*field_it)->getNbOfCoeffs();
          if(field_rank != dofPtr->getNbOfCoeffs()) {
            SETERRQ(PETSC_COMM_SELF,1,"data inconsistency");
          }
          double *tag_value;
          int tag_size;
          rval = mField.get_moab().tag_get_by_ptr(
            th,&node,1,(const void **)&tag_value,&tag_size
          ); CHKERRQ_MOAB(rval);
          if(tag_size != dofPtr->getNbOfCoeffs()) {
            SETERRQ(PETSC_COMM_SELF,1,"data inconsistency");
          }
          tag_value[dofPtr->getDofCoeffIdx()] = dofPtr->getFieldData();
        }
      }
      PetscFunctionReturn(0);
    }
    if(dofPtr->getEntType() != MBEDGE) {
      PetscFunctionReturn(0);
    }
    EntityHandle edge = dofPtr->getEnt();
    if(mField.get_moab().type_from_handle(edge)!=MBEDGE) {
      SETERRQ(PETSC_COMM_SELF,1,"this method works only elements which are type of MBEDGE");
    }

    int num_nodes;
    const EntityHandle* conn;
    rval = mField.get_moab().get_connectivity(edge,conn,num_nodes,false); CHKERRQ_MOAB(rval);
    if( (num_nodes != 2) && (num_nodes != 3) ) {
      SETERRQ(PETSC_COMM_SELF,1,"this method works only 4 node and 10 node tets");
    }
    if(num_nodes == 2) {
      PetscFunctionReturn(0);
    }

    if(dofPtr->getDofOrder()>=maxApproximationOrder) {
      SETERRQ(PETSC_COMM_SELF,1,"too big approximation order, increase constant max_ApproximationOrder");
    }
    double approx_val = 0;
    FieldApproximationBase base = dofPtr->getApproxBase();
    switch (base) {
      case AINSWORTH_COLE_BASE:
      approx_val = 0.25*L[dofPtr->getDofOrder()-2]*dofPtr->getFieldData();
      break;
      case LOBATTO_BASE:
      approx_val = 0.25*K[dofPtr->getDofOrder()-2]*dofPtr->getFieldData();
      break;
      default:
      SETERRQ(PETSC_COMM_SELF,MOFEM_NOT_IMPLEMENTED,"not yet implemented");
    }

    if(onCoords) {
      coords.resize(num_nodes*3);
      rval = mField.get_moab().get_coords(conn,num_nodes,&*coords.data().begin()); CHKERR_MOAB(rval);
      if(dofPtr->getEntDofIdx() == dofPtr->getDofCoeffIdx()) {
        //add only one when higher order terms present
        double ave_mid = (coords[3*0+dofPtr->getDofCoeffIdx()] + coords[3*1+dofPtr->getDofCoeffIdx()])*0.5;
        coords[2*3+dofPtr->getDofCoeffIdx()] = ave_mid;
      }
      coords[2*3+dofPtr->getDofCoeffIdx()] += approx_val;
      rval = mField.get_moab().set_coords(&conn[2],1,&coords[3*2]);  CHKERR_MOAB(rval);
    } else {
      int tag_size;
      double *tag_value[num_nodes];
      rval = mField.get_moab().tag_get_by_ptr(
        th,conn,num_nodes,(const void **)tag_value,&tag_size
      ); CHKERRQ_MOAB(rval);
      if(tag_size != dofPtr->getNbOfCoeffs()) {
        SETERRQ(PETSC_COMM_SELF,1,"data inconsistency");
      }
      if(dofPtr->getEntDofIdx() == dofPtr->getDofCoeffIdx()) {
        //add only one when higher order terms present
        double ave_mid = (tag_value[0][dofPtr->getDofCoeffIdx()] + tag_value[1][dofPtr->getDofCoeffIdx()])*0.5;
        tag_value[2][dofPtr->getDofCoeffIdx()] = ave_mid;
      }
      tag_value[2][dofPtr->getDofCoeffIdx()] += approx_val;
    }
    PetscFunctionReturn(0);
  }


};


}

#endif // __PROJECTION10NODECOORDSONFIELD_HPP__
