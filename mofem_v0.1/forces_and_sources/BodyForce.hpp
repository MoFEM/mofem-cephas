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

    PetscErrorCode doWork(
	int side,EntityType type,
	DataForcesAndSurcesCore::EntData &data) {
      PetscFunctionBegin;
      PetscErrorCode ierr;


      PetscFunctionReturn(0);
    }

  };

  opBodyForce opForce;
  int msId;

  BodyFroce(FieldInterface &_mField,Vec _F,const int _msId): 
    ForcesAndSurcesCore(_mField),opForce(),msId(_msId) {}

  const CubitMeshSets *cubit_meshset_ptr;
  EntityHandle meshset;
  PetscErrorCode getBodyElements() {
    PetscFunctionBegin;
    PetscErrorCode ierr;
    ierr = mField.get_Cubit_msId(msId,BlockSet|Block_BodyForcesSet,&cubit_meshset_ptr); CHKERRQ(ierr);
    //ierr = cubit_meshset_ptr->get_attribute_data_structure(opForce.blockData); CHKERRQ(ierr);
    meshset = cubit_meshset_ptr->get_meshset();
    PetscFunctionReturn(0);
  }

  PetscErrorCode operator()() {
    PetscFunctionBegin;
    PetscErrorCode ierr;
    ErrorCode rval;

    
    PetscFunctionReturn(0);
  }

};

}

#endif //__BODY_FORCE_HPP

