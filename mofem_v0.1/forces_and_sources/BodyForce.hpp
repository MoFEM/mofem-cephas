/* Copyright (C) 2013, Lukasz Kaczmarczyk (likask AT wp.pl)
 * --------------------------------------------------------------
 *
 * Test for linar elastic dynamics.
 *
 * This is not exactly procedure for linear elatic dynamics, since jacobian is
 * evaluated at every time step and snes procedure is involved. However it is
 * implemented like that, to test methodology for general nonlinear problem.
 *
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


#ifndef __BODY_FORCE_HPP
#define __BODY_FORCE_HPP

#include "ForcesAndSurcesCore.hpp"
extern "C" {
#include "gm_rule.h"
}

namespace MoFEM {

struct BodyFroce: public ForcesAndSurcesCore {

  struct opBodyForce: public DataOperator {

      Vec F;
      ublas::vector<double> &G_W;
      opBodyForce(Vec _F,ublas::vector<double> &_G_W):  F(_F),G_W(_G_W) {}

      Block_BodyForces blockData;

      double Volume;

      ublas::vector<FieldData> f;
      PetscErrorCode doWork(
	int side,EntityType type,
	DataForcesAndSurcesCore::EntData &data) {
      PetscFunctionBegin;
      PetscErrorCode ierr;

      f.resize(data.getIndices().size());
      bzero(&*f.data().begin(),data.getIndices().size()*sizeof(double));
      for(int gg = 0;gg<data.getN().size1();gg++) {
      }
      ierr = VecSetValues(F,data.getIndices().size(),&*data.getIndices().data().begin(),&*f.data().begin(),ADD_VALUES); CHKERRQ(ierr);

      PetscFunctionReturn(0);

      }
  };

  opBodyForce opForce;
  const int msId;

  ublas::vector<double> G_X,G_Y,G_Z,G_W;
  ublas::vector<double> coords;

  BodyFroce(FieldInterface &_mField,Vec _F,const int _msId): 
    ForcesAndSurcesCore(_mField),opForce(_F,G_W),msId(_msId) {}

  DataForcesAndSurcesCore data;
  DataForcesAndSurcesCore data_ho_gemetry;

  const CubitMeshSets *cubit_meshset_ptr;
  EntityHandle meshset;

  PetscErrorCode getBodyElements() {
    PetscFunctionBegin;
    PetscErrorCode ierr;
    ierr = mField.get_Cubit_msId(msId,BlockSet|Block_BodyForcesSet,&cubit_meshset_ptr); CHKERRQ(ierr);
    ierr = cubit_meshset_ptr->get_attribute_data_structure(opForce.blockData); CHKERRQ(ierr);
    meshset = cubit_meshset_ptr->get_meshset();
    PetscFunctionReturn(0);
  }

  PetscErrorCode operator()() {
    PetscFunctionBegin;
    PetscErrorCode ierr;
    ErrorCode rval;

    ierr = getEdgesSense(data); CHKERRQ(ierr);
    ierr = getFacesSense(data); CHKERRQ(ierr);

    ierr = getEdgesOrder(data); CHKERRQ(ierr);
    ierr = getFacesOrder(data); CHKERRQ(ierr);
    ierr = getOrderVolume(data); CHKERRQ(ierr);
    ierr = getRowNodesIndices(data,"DISPLACEMENT"); CHKERRQ(ierr);
    ierr = getEdgeRowIndices(data,"DISPLACEMENT"); CHKERRQ(ierr);
    ierr = getFacesRowIndices(data,"DISPLACEMENT"); CHKERRQ(ierr);
    ierr = getTetRowIndices(data,"DISPLACEMENT"); CHKERRQ(ierr);
    ierr = getFaceNodes(data); CHKERRQ(ierr);

    int order = 1;
    if(mField.check_field("MESH_NODE_POSITIONS")) {
      ierr = getNodesFieldData(data,"MESH_NODE_POSITIONS"); CHKERRQ(ierr);
      ierr = getEdgeFieldData(data,"MESH_NODE_POSITIONS"); CHKERRQ(ierr);
      ierr = getFacesFieldData(data,"MESH_NODE_POSITIONS"); CHKERRQ(ierr);
      ierr = getTetFieldData(data,"MESH_NODE_POSITIONS"); CHKERRQ(ierr);
    }


    int nb_gauss_pts = gm_rule_size(order,3);
    G_X.resize(nb_gauss_pts);
    G_Y.resize(nb_gauss_pts);
    G_Z.resize(nb_gauss_pts);
    G_W.resize(nb_gauss_pts);
    ierr = Grundmann_Moeller_integration_points_3D_TET(order,
	&*G_X.data().begin(),&*G_Y.data().begin(),&*G_Z.data().begin(),
	&*G_W.data().begin()); CHKERRQ(ierr);
    ierr = shapeTETFunctions_H1(data,
      &*G_X.data().begin(),
      &*G_X.data().begin(),
      &*G_X.data().begin(),nb_gauss_pts); CHKERRQ(ierr);

    EntityHandle ent = fe_ptr->get_ent();
    int num_nodes;
    const EntityHandle* conn;
    rval = mField.get_moab().get_connectivity(ent,conn,num_nodes,true); CHKERR_PETSC(rval);
    coords.resize(num_nodes*3);
    rval = mField.get_moab().get_coords(conn,num_nodes,&*coords.data().begin()); CHKERR_PETSC(rval);
    opForce.Volume = Shape_intVolumeMBTET(&*data.nOdes.getDiffN().data().begin(),&*coords.data().begin()); 

    /*opForce.invJac.resize(3,3);
    ierr = ShapeJacMBTET(&*data.nOdes.diffN.data().begin(),&*_coords_.begin(),&*opForce.invJac.data().begin()); CHKERRQ(ierr);
    ierr = Shape_invJac(&*opForce.invJac.data().begin()); CHKERRQ(ierr);*/

    //ierr = op.opNH1NH1(data); CHKERRQ(ierr);
    
    PetscFunctionReturn(0);
  }

};

}

#endif //__BODY_FORCE_HPP

