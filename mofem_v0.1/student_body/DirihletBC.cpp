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

#include "DirihletBC.hpp"

using namespace boost::numeric;

namespace MoFEM {

PetscErrorCode DisplacementBCFEMethodPreAndPostProc::iNitalize() {
  PetscFunctionBegin;
  if(map_zero_rows.empty()) {
    ParallelComm* pcomm = ParallelComm::get_pcomm(&mField.get_moab(),MYPCOMM_INDEX);
    for(_IT_CUBITMESHSETS_BY_BCDATA_TYPE_FOR_LOOP_(mField,NODESET|DISPLACEMENTSET,it)) {
	displacement_cubit_bc_data mydata;
	ierr = it->get_cubit_bc_data_structure(mydata); CHKERRQ(ierr);
	for(int dim = 0;dim<3;dim++) {
	  Range ents;
	  ierr = it->get_Cubit_msId_entities_by_dimension(mField.get_moab(),dim,ents,true); CHKERRQ(ierr);
	  if(dim>1) {
          Range _edges;
          ierr = mField.get_moab().get_adjacencies(ents,1,false,_edges,Interface::UNION); CHKERRQ(ierr);
          ents.insert(_edges.begin(),_edges.end());
	  }
        if(dim>0) {
          Range _nodes;
          rval = mField.get_moab().get_connectivity(ents,_nodes,true); CHKERR_PETSC(rval);
          ents.insert(_nodes.begin(),_nodes.end());
        }
	  for(Range::iterator eit = ents.begin();eit!=ents.end();eit++) {
	    for(_IT_NUMEREDDOFMOFEMENTITY_ROW_BY_NAME_ENT_PART_FOR_LOOP_(problem_ptr,fieldName,*eit,pcomm->rank(),dof)) {
	      if(dof->get_ent_type() == MBVERTEX) {
		if(dof->get_dof_rank() == 0 && mydata.data.flag1) {
		  map_zero_rows[dof->get_petsc_gloabl_dof_idx()] = mydata.data.value1;
		}
		if(dof->get_dof_rank() == 1 && mydata.data.flag2) {
		  map_zero_rows[dof->get_petsc_gloabl_dof_idx()] = mydata.data.value2;
		}
		if(dof->get_dof_rank() == 2 && mydata.data.flag3) {
		  map_zero_rows[dof->get_petsc_gloabl_dof_idx()] = mydata.data.value3;
		}
	      } else {
		if(dof->get_dof_rank() == 0 && mydata.data.flag1) {
		  map_zero_rows[dof->get_petsc_gloabl_dof_idx()] = 0;
		}
		if(dof->get_dof_rank() == 1 && mydata.data.flag2) {
		  map_zero_rows[dof->get_petsc_gloabl_dof_idx()] = 0;
		}
		if(dof->get_dof_rank() == 2 && mydata.data.flag3) {
		  map_zero_rows[dof->get_petsc_gloabl_dof_idx()] = 0;
		}
	      }
	    }
	  }
	}
    }
    dofsIndices.resize(map_zero_rows.size());
    dofsValues.resize(map_zero_rows.size());
    int ii = 0;
    map<DofIdx,FieldData>::iterator mit = map_zero_rows.begin();
    for(;mit!=map_zero_rows.end();mit++,ii++) { 
      dofsIndices[ii] = mit->first;
      dofsValues[ii] = mit->second;
    }
  }
  PetscFunctionReturn(0);
}

PetscErrorCode DisplacementBCFEMethodPreAndPostProc::preProcess() {
  PetscFunctionBegin;
  ierr = iNitalize(); CHKERRQ(ierr);

  if(snes_ctx == ctx_SNESNone && ts_ctx == ctx_TSNone) {
    ierr = VecSetValues(snes_x,dofsIndices.size(),&dofsIndices[0],&dofsValues[0],INSERT_VALUES); CHKERRQ(ierr);
    ierr = VecAssemblyBegin(snes_x); CHKERRQ(ierr);
    ierr = VecAssemblyEnd(snes_x); CHKERRQ(ierr);
  }

  switch(snes_ctx) {
    case ctx_SNESNone: {} 
    break;
    case ctx_SNESSetFunction: {
      ierr = VecSetValues(snes_x,dofsIndices.size(),&dofsIndices[0],&dofsValues[0],INSERT_VALUES); CHKERRQ(ierr);
      ierr = VecAssemblyBegin(snes_x); CHKERRQ(ierr);
      ierr = VecAssemblyEnd(snes_x); CHKERRQ(ierr);
    }
    break;
    case ctx_SNESSetJacobian: {
    }
    break;
    default:
	SETERRQ(PETSC_COMM_SELF,1,"unknown snes stage");
  }

  switch(ts_ctx) {
    case ctx_TSNone: {}
    break;
    case ctx_TSSetIFunction: {
      ierr = VecSetValues(ts_u,dofsIndices.size(),&dofsIndices[0],&dofsValues[0],INSERT_VALUES); CHKERRQ(ierr);
      ierr = VecAssemblyBegin(ts_u); CHKERRQ(ierr);
      ierr = VecAssemblyEnd(ts_u); CHKERRQ(ierr);
    }
    break;
    case ctx_TSSetIJacobian: {
    }
    break;
    default:
	SETERRQ(PETSC_COMM_SELF,1,"unknown snes stage");
  }


  PetscFunctionReturn(0);
}

PetscErrorCode DisplacementBCFEMethodPreAndPostProc::postProcess() {
  PetscFunctionBegin;
  if(snes_ctx == ctx_SNESNone && ts_ctx == ctx_TSNone) {
    ierr = MatAssemblyBegin(*snes_B,MAT_FINAL_ASSEMBLY); CHKERRQ(ierr);
    ierr = MatAssemblyEnd(*snes_B,MAT_FINAL_ASSEMBLY); CHKERRQ(ierr);
    ierr = MatZeroRowsColumns(*snes_B,dofsIndices.size(),&dofsIndices[0],1,PETSC_NULL,PETSC_NULL); CHKERRQ(ierr);
    ierr = VecAssemblyBegin(snes_f); CHKERRQ(ierr);
    ierr = VecAssemblyEnd(snes_f); CHKERRQ(ierr);
    for(vector<int>::iterator vit = dofsIndices.begin();vit!=dofsIndices.end();vit++) {
      ierr = VecSetValue(snes_f,*vit,0,INSERT_VALUES); CHKERRQ(ierr);
    }
    ierr = VecAssemblyBegin(snes_f); CHKERRQ(ierr);
    ierr = VecAssemblyEnd(snes_f); CHKERRQ(ierr);
  }

  switch(snes_ctx) {
    case ctx_SNESNone: {}
    break;
    case ctx_SNESSetFunction: {
      ierr = VecAssemblyBegin(snes_f); CHKERRQ(ierr);
      ierr = VecAssemblyEnd(snes_f); CHKERRQ(ierr);
      for(vector<int>::iterator vit = dofsIndices.begin();vit!=dofsIndices.end();vit++) {
	ierr = VecSetValue(snes_f,*vit,0,INSERT_VALUES); CHKERRQ(ierr);
      }
      ierr = VecAssemblyBegin(snes_f); CHKERRQ(ierr);
      ierr = VecAssemblyEnd(snes_f); CHKERRQ(ierr);
    }
    break;
    case ctx_SNESSetJacobian: {
      ierr = MatAssemblyBegin(*snes_B,MAT_FINAL_ASSEMBLY); CHKERRQ(ierr);
      ierr = MatAssemblyEnd(*snes_B,MAT_FINAL_ASSEMBLY); CHKERRQ(ierr);
      ierr = MatZeroRowsColumns(*snes_B,dofsIndices.size(),&dofsIndices[0],1,PETSC_NULL,PETSC_NULL); CHKERRQ(ierr);
    }
    break;
    default:
	SETERRQ(PETSC_COMM_SELF,1,"unknown snes stage");
  }

  switch(ts_ctx) {
    case ctx_TSNone: {}
    break;
    case ctx_TSSetIFunction: {
      ierr = VecAssemblyBegin(ts_F); CHKERRQ(ierr);
      ierr = VecAssemblyEnd(ts_F); CHKERRQ(ierr);
      for(vector<int>::iterator vit = dofsIndices.begin();vit!=dofsIndices.end();vit++) {
	ierr = VecSetValue(ts_F,*vit,0,INSERT_VALUES); CHKERRQ(ierr);
      }
      ierr = VecAssemblyBegin(ts_F); CHKERRQ(ierr);
      ierr = VecAssemblyEnd(ts_F); CHKERRQ(ierr);
    }
    break;
    case ctx_TSSetIJacobian: {
      ierr = MatAssemblyBegin(*ts_B,MAT_FINAL_ASSEMBLY); CHKERRQ(ierr);
      ierr = MatAssemblyEnd(*ts_B,MAT_FINAL_ASSEMBLY); CHKERRQ(ierr);
      ierr = MatZeroRowsColumns(*ts_B,dofsIndices.size(),&dofsIndices[0],1,PETSC_NULL,PETSC_NULL); CHKERRQ(ierr);
    }
    break;
    default:
	SETERRQ(PETSC_COMM_SELF,1,"unknown snes stage");
  }

  PetscFunctionReturn(0);
}

PetscErrorCode SpatialPositionsBCFEMethodPreAndPostProc::iNitalize() {
  PetscFunctionBegin;
  if(map_zero_rows.empty()) {
    ParallelComm* pcomm = ParallelComm::get_pcomm(&mField.get_moab(),MYPCOMM_INDEX);
    for(_IT_CUBITMESHSETS_BY_BCDATA_TYPE_FOR_LOOP_(mField,NODESET|DISPLACEMENTSET,it)) {
	displacement_cubit_bc_data mydata;
	ierr = it->get_cubit_bc_data_structure(mydata); CHKERRQ(ierr);
	for(int dim = 0;dim<3;dim++) {
	  Range ents;
	  ierr = it->get_Cubit_msId_entities_by_dimension(mField.get_moab(),dim,ents,true); CHKERRQ(ierr);
	  if(dim>1) {
          Range _edges;
          ierr = mField.get_moab().get_adjacencies(ents,1,false,_edges,Interface::UNION); CHKERRQ(ierr);
          ents.insert(_edges.begin(),_edges.end());
	  }
        if(dim>0) {
          Range _nodes;
          rval = mField.get_moab().get_connectivity(ents,_nodes,true); CHKERR_PETSC(rval);
          ents.insert(_nodes.begin(),_nodes.end());
        }
	  for(Range::iterator eit = ents.begin();eit!=ents.end();eit++) {
	    for(_IT_NUMEREDDOFMOFEMENTITY_ROW_BY_NAME_ENT_PART_FOR_LOOP_(problem_ptr,fieldName,*eit,pcomm->rank(),dof)) {
	      if(dof->get_ent_type() == MBVERTEX) {
		EntityHandle node = dof->get_ent();
		cOords.resize(3);
		rval = mField.get_moab().get_coords(&node,1,&*cOords.data().begin()); CHKERR_PETSC(rval);
		if(dof->get_dof_rank() == 0 && mydata.data.flag1) {
		  map_zero_rows[dof->get_petsc_gloabl_dof_idx()] = cOords[0]+mydata.data.value1;
		}
		if(dof->get_dof_rank() == 1 && mydata.data.flag2) {
		  map_zero_rows[dof->get_petsc_gloabl_dof_idx()] = cOords[1]+mydata.data.value2;
		}
		if(dof->get_dof_rank() == 2 && mydata.data.flag3) {
		  map_zero_rows[dof->get_petsc_gloabl_dof_idx()] = cOords[2]+mydata.data.value3;
		}
	      } else {
		if(dof->get_dof_rank() == 0 && mydata.data.flag1) {
		  map_zero_rows[dof->get_petsc_gloabl_dof_idx()] = 0;
		}
		if(dof->get_dof_rank() == 1 && mydata.data.flag2) {
		  map_zero_rows[dof->get_petsc_gloabl_dof_idx()] = 0;
		}
		if(dof->get_dof_rank() == 2 && mydata.data.flag3) {
		  map_zero_rows[dof->get_petsc_gloabl_dof_idx()] = 0;
		}
	      }
	    }
	    for(vector<string>::iterator fit = fixFields.begin();fit!=fixFields.end();fit++) {
	      for(_IT_NUMEREDDOFMOFEMENTITY_ROW_BY_NAME_ENT_PART_FOR_LOOP_(problem_ptr,*fit,*eit,pcomm->rank(),dof)) {
		map_zero_rows[dof->get_petsc_gloabl_dof_idx()] = dof->get_FieldData();
	      }
	    }
	  }
	}
    }
    dofsIndices.resize(map_zero_rows.size());
    dofsValues.resize(map_zero_rows.size());
    int ii = 0;
    map<DofIdx,FieldData>::iterator mit = map_zero_rows.begin();
    for(;mit!=map_zero_rows.end();mit++,ii++) { 
      dofsIndices[ii] = mit->first;
      dofsValues[ii] = mit->second;
    }
  }
  PetscFunctionReturn(0);
}

PetscErrorCode TemperatureBCFEMethodPreAndPostProc::iNitalize() {
  PetscFunctionBegin;
  if(map_zero_rows.empty()) {
    ParallelComm* pcomm = ParallelComm::get_pcomm(&mField.get_moab(),MYPCOMM_INDEX);

    for(_IT_CUBITMESHSETS_BY_BCDATA_TYPE_FOR_LOOP_(mField,NODESET|TEMPERATURESET,it)) {
      temperature_cubit_bc_data mydata;
      ierr = it->get_cubit_bc_data_structure(mydata); CHKERRQ(ierr);
      for(int dim = 0;dim<3;dim++) {
        Range ents;
        ierr = it->get_Cubit_msId_entities_by_dimension(mField.get_moab(),dim,ents,true); CHKERRQ(ierr);
        if(dim>1) {
  	Range _edges;
  	ierr = mField.get_moab().get_adjacencies(ents,1,false,_edges,Interface::UNION); CHKERRQ(ierr);
  	ents.insert(_edges.begin(),_edges.end());
        }
        if(dim>0) {
	  Range _nodes;
	  rval = mField.get_moab().get_connectivity(ents,_nodes,true); CHKERR_PETSC(rval);
	  ents.insert(_nodes.begin(),_nodes.end());
        }
        for(Range::iterator eit = ents.begin();eit!=ents.end();eit++) {
  	for(_IT_NUMEREDDOFMOFEMENTITY_ROW_BY_NAME_ENT_PART_FOR_LOOP_(problem_ptr,fieldName,*eit,pcomm->rank(),dof)) {
  	  if(dof->get_ent_type() == MBVERTEX) {
  	    map_zero_rows[dof->get_petsc_gloabl_dof_idx()] = mydata.data.value1;
  	  } else {
  	    map_zero_rows[dof->get_petsc_gloabl_dof_idx()] = 0;
  	  }
  	}
        }
      }
    }
    dofsIndices.resize(map_zero_rows.size());
    dofsValues.resize(map_zero_rows.size());
    int ii = 0;
    map<DofIdx,FieldData>::iterator mit = map_zero_rows.begin();
    for(;mit!=map_zero_rows.end();mit++,ii++) { 
      dofsIndices[ii] = mit->first;
      dofsValues[ii] = mit->second;
    }

  }
  PetscFunctionReturn(0);
}

PetscErrorCode FixMaterialPoints::iNitalize() {
  PetscFunctionBegin;
  ParallelComm* pcomm = ParallelComm::get_pcomm(&mField.get_moab(),MYPCOMM_INDEX);
  if(map_zero_rows.empty()) {
    for(vector<string>::iterator fit = fieldNames.begin();fit!=fieldNames.end();fit++) {
      for(Range::iterator eit = eNts.begin();eit!=eNts.end();eit++) {
	for(_IT_NUMEREDDOFMOFEMENTITY_ROW_BY_NAME_ENT_PART_FOR_LOOP_(problem_ptr,*fit,*eit,pcomm->rank(),dof)) {
	 map_zero_rows[dof->get_petsc_gloabl_dof_idx()] = 0;
	}
      }
    }
    dofsIndices.resize(map_zero_rows.size());
    dofsValues.resize(map_zero_rows.size());
    int ii = 0;
    map<DofIdx,FieldData>::iterator mit = map_zero_rows.begin();
    for(;mit!=map_zero_rows.end();mit++,ii++) { 
      dofsIndices[ii] = mit->first;
      dofsValues[ii] = mit->second;
    }
  }
  PetscFunctionReturn(0);
}

PetscErrorCode FixMaterialPoints::preProcess() {
    PetscFunctionBegin;
    ierr = iNitalize(); CHKERRQ(ierr);
    PetscFunctionReturn(0);
  }

}

