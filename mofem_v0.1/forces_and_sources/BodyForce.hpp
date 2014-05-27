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

  const int msId;
  CoreForcesAndSurces(FiedInterface& _mField,const int _msId):
    CoreForcesAndSurces(mField),msId(_msId) {};

  struct opBodyForce: public dataOperator {

      Block_BodyForces data;

      ublas::vector<FieldData> f;
      PetscErrorCode doWork(
	int side,int side,
	EntityType type
	ublas::vector<DofIdx> &indices
	ublas::matrix<FieldData> &N) {

      PetscFunctionBegin;

      f.resize(indices.size());
      for(int gg = 0;gg<N.size1();gg++) {
	
      }


      PetscFunctionReturn(0);

      }
  }

  opBodyForce op;

  const CubitMeshSets *cubit_meshset_ptr;
  EntityHandle meshset;

  PetscErroCode getBodyElements() {
    PetscFunctionBegin;
    ierr = mField.get_Cubit_msId(msId,BlockSet|Body_Force,&cubit_meshset_ptr); CHKERRQ(ierr);
    ierr = cubit_meshset_ptr->get_attribute_data_structure(op.data); CHKERRQ(ierr);
    meshset = cubit_meshset_ptr->get_meshset();
    PetscFunctionReturn(0);
  }

  dataForcesAndSurcesCore data;

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

    ierr = shapeTETFunctions_H1(data,G_TET_X4,G_TET_Y4,G_TET_Z4,4); CHKERRQ(ierr);

    ierr = op.opNH1NH1(data); CHKERRQ(ierr);

    PetscFunctionReturn(0);
  }

};

}

#endif //__BODY_FORCE_HPP

