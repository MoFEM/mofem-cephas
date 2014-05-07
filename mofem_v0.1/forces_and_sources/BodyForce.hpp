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

#include "CoreForcesAndSurces.hpp"

namespace MoFEM {

struct BodyFroce: public CoreForcesAndSurces {

  struct opBodyForce: public dataOperator {

      ublas::vector<double> &G_W;
      opBodyForce(Vec _F,blas::vector<double> &_G_W): F(_F),G_W(_G_W) {}

      double Volume;
      ublas::matrix<FieldData> invJac;
      ublas::matrix<FieldData> diffNinvJac;

      ublas::vector<FieldData> f;
      PetscErrorCode doWork(
	int side,EntityType type,
	dataForcesAndSurcesCore::entData &data) {

      PetscFunctionBegin;

      f.resize(indices.size());
      bzero(&*f.data().begin(),indices.size()*sizeof(double));
      for(int gg = 0;gg<N.size1();gg++) {
	double val = G_W[gg]*Volume*data.density;
      }
      ierr = VecSetValues(F,indices.size(),&indices.data().begin(),&*f.data().begin(),ADD_VALUES); CHKERRQ(ierr);

      PetscFunctionReturn(0);

      }
  }

  struct opJacobian: public dataOperator {

      ublas::matrix<FieldData> &invJac;
      ublas::matrix<FieldData> Jac;

      opBodyForce(ublas::matrix<FieldData>  _invJac): invJac(_invJac) {}

      ublas::vector<FieldData> f;
      PetscErrorCode doWork(
	int side,EntityType type,
	dataForcesAndSurcesCore::entData &data) {
      PetscFunctionBegin;

      invJac.resize(data.N.size1(),9);

      Jac.resize(data.N.size1(),9);
      if(type == MBVERTEX) {
	bzero(&*Jac.data().begin(),Jac.size1()*Jac.size2()*sizeof(double));
      }
      typedef ublas::array_adaptor<FieldData> storage_t;

      int nb_dofs = data.N.size2();
      for(int gg = 0;gg<N.size1();gg++) {
	storage_t st_jac(3*3,&Jac(gg,0));
  	ublas::matrix<FieldData,ublas::row_major,storage_t> jac_at_GaussPt(3,3,st);
	for(int i = 0;i<3;i++) {
	  for(int j = 0;j<3;j++) {
	    jac_at_GaussPt(i,j) +=
	      cblas_ddot(nb_dofs,&data.diffN(gg,j),3,&data.fieldData(i),3);
	  }
	}
      }

      PetscFunctionReturn(0);
      }

  }

  opBodyForce opForce;
  opJacobiab opJac;

  ublas::vector<double> G_X,G_Y,G_Z,G_W;
  ublas::vector<double> coords;

  BodyFroce(mField &_mField,Vec _F): 
    mField(_mField),opForce(F,G_W),
    opJac(opForce.invJac) {}

  dataForcesAndSurcesCore data;
  dataForcesAndSurcesCore data_ho_gemetry;

  const CubitMeshSets *cubit_meshset_ptr;
  EntityHandle meshset;

  PetscErroCode getBodyElements() {
    PetscFunctionBegin;
    ierr = mField.get_Cubit_msId(msId,BlockSet|Body_Force,&cubit_meshset_ptr); CHKERRQ(ierr);
    ierr = cubit_meshset_ptr->get_attribute_data_structure(op.data); CHKERRQ(ierr);
    meshset = cubit_meshset_ptr->get_meshset();
    PetscFunctionReturn(0);
  }

  PetscErrorCode operator()() {
    PetscFunctionBegin;

    ierr = getEdgesSense(data); CHKERRQ(ierr);
    ierr = getFacesSense(data); CHKERRQ(ierr);

    ierr = getEdgesOrder(data,"DISPLACEMENT"); CHKERRQ(ierr);
    ierr = getFacesOrder(data,"DISPLACEMENT"); CHKERRQ(ierr);
    ierr = getOrderVolume(data,"DISPLACEMENT"); CHKERRQ(ierr);
    ierr = getRowNodesIndices(data,"DISPLACEMENT"); CHKERRQ(ierr);
    ierr = getEdgeRowIndices(data,"DISPLACEMENT"); CHKERRQ(ierr);
    ierr = getFacesRowIndices(data,"DISPLACEMENT"); CHKERRQ(ierr);
    ierr = getTetRowIndices(data,"DISPLACEMENT"); CHKERRQ(ierr);
    ierr = getFaceNodes(data); CHKERRQ(ierr);

    int order = 1;
    order = max(order,*max_element(data_.edgesOrder.begin(),data.edgesOrder.end()));
    order = max(order,*max_element(data.facesOrder.begin(),data.facesOrder.end()));
    order = max(order,data.volumeOrder);
  
    int ho_geom_order = 1;
    if(check_field("MESH_NODE_POSITIONS")) {
    
      ierr = getEdgesOrder(data_ho_gemetry,"MESH_NODE_POSITIONS"); CHKERRQ(ierr);
      ierr = getFacesOrder(data_ho_gemetry,"MESH_NODE_POSITIONS"); CHKERRQ(ierr);
      ierr = getOrderVolume(data_ho_gemetry,"MESH_NODE_POSITIONS"); CHKERRQ(ierr);

      ho_geom_order = max(ho_geom_order,*max_element(data_ho_gemetry_.edgesOrder.begin(),data_ho_gemetry.edgesOrder.end()));
      ho_geom_order = max(ho_geom_order,*max_element(data_ho_gemetry.facesOrder.begin(),data_ho_gemetry.facesOrder.end()));
      ho_geom_order = max(ho_geom_order,data_ho_gemetry.volumeOrder);

      ierr = getNodesFieldData(data,"MESH_NODE_POSITIONS"); CHKERRQ(ierr);
      ierr = getEdgeFieldData(data,"MESH_NODE_POSITIONS"); CHKERRQ(ierr);
      ierr = getFacesFieldData(data,"MESH_NODE_POSITIONS"); CHKERRQ(ierr);
      ierr = getTetFieldData(data,"MESH_NODE_POSITIONS"); CHKERRQ(ierr);
    }

    order = max(order,ho_geom_order);

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
    rval = moab.get_connectivity(eny,conn,num_nodes,true); CHKERR_PETSC(rval);
    coords.resize(num_nodes*3);
    rval = moab.get_coords(conn,num_nodes,&*coords.data().begin()); CHKERR_PETSC(rval);
    opForce.Volume = Shape_intVolumeMBTET(diffNTET,&*coords.data().begin()); 

    opForce.invJac.resize(3,3);
    ierr = ShapeJacMBTET(diffNTET,&*_coords_.begin(),&*opForce.invJac.data().begin()); CHKERRQ(ierr);
    ierr = Shape_invJac(&*opForce.invJac.data().begin()); CHKERRQ(ierr);
    opForce.diffNinvJac.resize(4,3);
    ierr = ShapeDiffMBTETinvJ(
      &*data.nOdes.diffN.data().begin(),
      &*opForce.invJac.data().begin(),&*opForce.diffNinvJac.data().begin()); CHKERRQ(ierr);


    ierr = op.opNH1NH1(data); CHKERRQ(ierr);
    
    PetscFunctionReturn(0);
  }

};

}

#endif //__BODY_FORCE_HPP

