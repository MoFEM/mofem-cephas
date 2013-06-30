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

#include "moabFEMethod_ComplexForLazy.hpp"
#include "FEM.h"
#include "complex_for_lazy.h"

namespace MoFEM {

PetscErrorCode FEMethod_ComplexForLazy::GetIndices() {
  PetscFunctionBegin;
  switch (fe_ent_ptr->get_ent_type()) {
    case MBTET: {
      RowGlob.resize(1+6+4+1);
      ColGlob.resize(1+6+4+1);
      rowNMatrices.resize(1+6+4+1);
      rowDiffNMatrices.resize(1+6+4+1);
      rowBMatrices.resize(1+6+4+1);
      colNMatrices.resize(1+6+4+1);
      colDiffNMatrices.resize(1+6+4+1);
      colBMatrices.resize(1+6+4+1);
      try {
      //nodes
      ierr = GetRowIndices("SPATIAL_POSITION",RowGlob[0]); CHKERRQ(ierr);
      ierr = GetColIndices("SPATIAL_POSITION",ColGlob[0]); CHKERRQ(ierr);
      //row
      ierr = GetGaussRowNMatrix("SPATIAL_POSITION",rowNMatrices[0]); CHKERRQ(ierr);
      ierr = GetGaussRowDiffNMatrix("SPATIAL_POSITION",rowDiffNMatrices[0]); CHKERRQ(ierr);
      ierr = MakeBMatrix3D("SPATIAL_POSITION",rowDiffNMatrices[0],rowBMatrices[0]);  CHKERRQ(ierr);
      //col
      ierr = GetGaussColNMatrix("SPATIAL_POSITION",colNMatrices[0]); CHKERRQ(ierr);
      ierr = GetGaussColDiffNMatrix("SPATIAL_POSITION",colDiffNMatrices[0]); CHKERRQ(ierr);
      ierr = MakeBMatrix3D("SPATIAL_POSITION",colDiffNMatrices[0],colBMatrices[0]);  CHKERRQ(ierr);
      //edges
      int ee = 0;
      for(;ee<6;ee++) { //edges matrices
	ierr = GetRowIndices("SPATIAL_POSITION",MBEDGE,RowGlob[1+ee],ee); CHKERRQ(ierr);
	if(RowGlob[1+ee].size()!=0) {
	  ierr = GetGaussRowNMatrix("SPATIAL_POSITION",MBEDGE,rowNMatrices[1+ee],ee); CHKERRQ(ierr);
	  ierr = GetGaussRowDiffNMatrix("SPATIAL_POSITION",MBEDGE,rowDiffNMatrices[1+ee],ee); CHKERRQ(ierr);
	  ierr = MakeBMatrix3D("SPATIAL_POSITION",rowDiffNMatrices[1+ee],rowBMatrices[1+ee]);  CHKERRQ(ierr);
	}
	ierr = GetColIndices("SPATIAL_POSITION",MBEDGE,ColGlob[1+ee],ee); CHKERRQ(ierr);
	if(ColGlob[1+ee].size()>0) {
	  ierr = GetGaussColNMatrix("SPATIAL_POSITION",MBEDGE,colNMatrices[1+ee],ee); CHKERRQ(ierr);
	  ierr = GetGaussColDiffNMatrix("SPATIAL_POSITION",MBEDGE,colDiffNMatrices[1+ee],ee); CHKERRQ(ierr);
	  ierr = MakeBMatrix3D("SPATIAL_POSITION",colDiffNMatrices[1+ee],colBMatrices[1+ee]);  CHKERRQ(ierr);
	}
      }
      assert(ee == 6);
      //faces
      int ff = 0;
      for(;ff<4;ff++) { //faces matrices
	ierr = GetRowIndices("SPATIAL_POSITION",MBTRI,RowGlob[1+ee+ff],ff); CHKERRQ(ierr);
	if(RowGlob[1+ee+ff].size()!=0) {
	  ierr = GetGaussRowNMatrix("SPATIAL_POSITION",MBTRI,rowNMatrices[1+ee+ff],ff); CHKERRQ(ierr);
	  ierr = GetGaussRowDiffNMatrix("SPATIAL_POSITION",MBTRI,rowDiffNMatrices[1+ee+ff],ff); CHKERRQ(ierr);
	  ierr = MakeBMatrix3D("SPATIAL_POSITION",rowDiffNMatrices[1+ee+ff],rowBMatrices[1+ee+ff]);  CHKERRQ(ierr);
	}
	ierr = GetColIndices("SPATIAL_POSITION",MBTRI,ColGlob[1+ee+ff],ff); CHKERRQ(ierr);
	if(ColGlob[1+ee+ff].size()!=0) {
	  ierr = GetGaussColNMatrix("SPATIAL_POSITION",MBTRI,colNMatrices[1+ee+ff],ff); CHKERRQ(ierr);
	  ierr = GetGaussColDiffNMatrix("SPATIAL_POSITION",MBTRI,colDiffNMatrices[1+ee+ff],ff); CHKERRQ(ierr);
	  ierr = MakeBMatrix3D("SPATIAL_POSITION",colDiffNMatrices[1+ee+ff],colBMatrices[1+ee+ff]);  CHKERRQ(ierr);
	}
      }
      assert(ff == 4);
      //volumes
      ierr = GetRowIndices("SPATIAL_POSITION",MBTET,RowGlob[1+ee+ff]); CHKERRQ(ierr);
      if(RowGlob[1+ee+ff].size()!=0) {
	ierr = GetGaussRowNMatrix("SPATIAL_POSITION",MBTET,rowNMatrices[1+ee+ff]); CHKERRQ(ierr);
	ierr = GetGaussRowDiffNMatrix("SPATIAL_POSITION",MBTET,rowDiffNMatrices[1+ee+ff]); CHKERRQ(ierr);
	ierr = MakeBMatrix3D("SPATIAL_POSITION",rowDiffNMatrices[1+ee+ff],rowBMatrices[1+ee+ff]);  CHKERRQ(ierr);
      }
      ierr = GetColIndices("SPATIAL_POSITION",MBTET,ColGlob[1+ee+ff]); CHKERRQ(ierr);
      if(ColGlob[1+ee+ff].size()!=0) {
	ierr = GetGaussColNMatrix("SPATIAL_POSITION",MBTET,colNMatrices[1+ee+ff]); CHKERRQ(ierr);
	ierr = GetGaussColDiffNMatrix("SPATIAL_POSITION",MBTET,colDiffNMatrices[1+ee+ff]); CHKERRQ(ierr);
	ierr = MakeBMatrix3D("SPATIAL_POSITION",colDiffNMatrices[1+ee+ff],colBMatrices[1+ee+ff]);  CHKERRQ(ierr);
      }
      } catch (const char* msg) {
	SETERRQ(PETSC_COMM_SELF,1,msg);
      } 
    }
    break;
    default:
    SETERRQ(PETSC_COMM_SELF,1,"no implemented");
  }
  PetscFunctionReturn(0);
}
PetscErrorCode FEMethod_ComplexForLazy::GetTangent() {
  PetscFunctionBegin;

  PetscFunctionReturn(0);
}
PetscErrorCode FEMethod_ComplexForLazy::GerResidual() {
  PetscFunctionBegin;
  PetscFunctionReturn(0);
}

}

