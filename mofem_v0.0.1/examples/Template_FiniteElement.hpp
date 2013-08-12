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

#include "moabField.hpp"
#include "moabField_Core.hpp"
#include "moabFEMethod_UpLevelStudent.hpp"
#include "cholesky.hpp"
#include <petscksp.h>

#include "ElasticFEMethod.hpp"
#include "PostProcVertexMethod.hpp"
#include "PostProcDisplacementAndStrainOnRefindedMesh.hpp"
#include "ElasticFEMethodForInterface.hpp"

#ifdef __cplusplus
extern "C" {
#endif
#include <petsc-private/snesimpl.h>
#ifdef __cplusplus
}
#endif

namespace MoFEM {

struct TemplateFEMethodSNES: public FEMethod_UpLevelStudent {

  TemplateFEMethodSNES(Interface& _moab): FEMethod_UpLevelStudent(_moab) {};
  TemplateFEMethodSNES(Interface& _moab,BaseDirihletBC *_dirihlet_bc_method_ptr,int verbose = 0): FEMethod_UpLevelStudent(_moab,_dirihlet_bc_method_ptr,verbose) {}

  PetscErrorCode preProcess() {
    PetscFunctionBegin;

    g_NTET.resize(4*45);
    ShapeMBTET(&g_NTET[0],G_TET_X45,G_TET_Y45,G_TET_Z45,45);
    G_W_TET = G_TET_W45;
    g_NTRI.resize(3*13);
    ShapeMBTRI(&g_NTRI[0],G_TRI_X13,G_TRI_Y13,13); 
    G_W_TRI = G_TRI_W13;

    switch(snes_ctx) {
      case ctx_SNESNone: {
      }
      break;
      case ctx_SNESSetFunction: { 
      }
      break;
      case ctx_SNESSetJacobian: {
      }
      break;
      default:
	SETERRQ(PETSC_COMM_SELF,1,"not implemented");
    }

    PetscFunctionReturn(0);
  }


  PetscErrorCode operator()() {
      PetscFunctionBegin;
      ierr = OpStudentStart_TET(g_NTET); CHKERRQ(ierr);

      switch(snes_ctx) {
	case ctx_SNESNone: {
	}
	break;
	case ctx_SNESSetJacobian: {

	  for(int rr = 0;rr<row_mat;rr++) {
	    if(RowGlob[rr].size()==0) continue;
	    ierr = VecSetValues(snes_f,/*LOCAL_VECTOR_TO_ASSEMBLE*/,ADD_VALUES); CHKERRQ(ierr);
	  }

	break;
	case ctx_SNESSetFunction: { 

	  for(int rr = 0;rr<row_mat;rr++) {
	    if(RowGlob[rr].size()==0) continue;
	    for(int cc = 0;cc<col_mat;cc++) {
	      if(ColGlob[cc].size()==0) continue;
		// ****
		ierr = MatSetValues(*snes_B,/*LOCAL_MATRIX_TO_ASSEMBLE*/,ADD_VALUES); CHKERRQ(ierr);
		// ****
	    }
	  }
	

	}
	break;
	default:
	  SETERRQ(PETSC_COMM_SELF,1,"not implemented");
      }

      ierr = OpStudentEnd(); CHKERRQ(ierr);
      PetscFunctionReturn(0);
  }


  PetscErrorCode postProcess() {
    PetscFunctionBegin;

    switch(snes_ctx) {
      case ctx_SNESNone: {
      }
      case ctx_SNESSetFunction: { 
      }
      break;
      case ctx_SNESSetJacobian: {

      }
      break;
      default:
	SETERRQ(PETSC_COMM_SELF,1,"not implemented");
    }

    PetscFunctionReturn(0);
  }

  private:

    vector<double> g_NTET,g_NTRI;
    const double* G_W_TET;
    const double* G_W_TRI;

};


}
