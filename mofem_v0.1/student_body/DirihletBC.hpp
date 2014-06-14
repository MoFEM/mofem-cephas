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

#ifndef __MOABFEMETHOD_DIRIHLETBC_HPP__
#define __MOABFEMETHOD_DIRIHLETBC_HPP__

#include "FieldInterface.hpp"

#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/matrix_proxy.hpp>
#include <boost/numeric/ublas/io.hpp>

using namespace boost::numeric;

namespace MoFEM {

struct DisplacementBCFEMethodPreAndPostProc: public FieldInterface::FEMethod {

  FieldInterface& mField;
  const string field_name;

  DisplacementBCFEMethodPreAndPostProc(
    FieldInterface& _mField,const string &_field_name,
    Mat &_Aij,Vec _X,Vec _F): mField(_mField),field_name(_field_name) {
    snes_B = &_Aij;
    snes_x = _X;
    snes_f = _F;
    ts_B = &_Aij;
  };

  PetscErrorCode ierr;
  ErrorCode rval;

  map<DofIdx,FieldData> map_zero_rows;
  vector<int> dofsIndices;
  vector<double> dofsValues;

  virtual PetscErrorCode iNitalize() {
    PetscFunctionBegin;
    if(map_zero_rows.empty()) {
      ParallelComm* pcomm = ParallelComm::get_pcomm(&mField.get_moab(),MYPCOMM_INDEX);
      for(_IT_CUBITMESHSETS_BY_BCDATA_TYPE_FOR_LOOP_(mField,NodeSet|DisplacementSet,it)) {
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
	    for(_IT_NUMEREDDOFMOFEMENTITY_COL_BY_ENT_FOR_LOOP_(problem_ptr,*eit,dof)) {
	      if(dof->get_name()!=field_name) continue;
	      if(dof->get_part()!=pcomm->rank()) continue;
	      if(dof->get_dof_rank() == 0 && mydata.data.flag1) {
		map_zero_rows[dof->get_petsc_gloabl_dof_idx()] = mydata.data.value1;
	      }
	      if(dof->get_dof_rank() == 1 && mydata.data.flag2) {
		map_zero_rows[dof->get_petsc_gloabl_dof_idx()] = mydata.data.value2;
	      }
	      if(dof->get_dof_rank() == 2 && mydata.data.flag3) {
		map_zero_rows[dof->get_petsc_gloabl_dof_idx()] = mydata.data.value3;
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

  PetscErrorCode preProcess() {
    PetscFunctionBegin;
    ierr = iNitalize(); CHKERRQ(ierr);
    ierr = VecSetValues(snes_x,dofsIndices.size(),&dofsIndices[0],&dofsValues[0],INSERT_VALUES); CHKERRQ(ierr);
    ierr = VecAssemblyBegin(snes_x); CHKERRQ(ierr);
    ierr = VecAssemblyEnd(snes_x); CHKERRQ(ierr);
    PetscFunctionReturn(0);
  }

  PetscErrorCode postProcess() {
    PetscFunctionBegin;
    switch(snes_ctx) {
      case ctx_SNESNone: {
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
    PetscFunctionReturn(0);
  }

};


/** 
 * \brief The student user interface for Dirihlet boundary conditions
 * 
*/
struct BaseDirihletBC {

  BaseDirihletBC();


  virtual PetscErrorCode SetDirihletBC_to_ElementIndiciesRow(
    FieldInterface::FEMethod *fe_method_ptr,vector<vector<DofIdx> > &RowGlobDofs,vector<DofIdx>& DirihletBC);
  virtual PetscErrorCode SetDirihletBC_to_ElementIndiciesCol(
    FieldInterface::FEMethod *fe_method_ptr,vector<vector<DofIdx> > &ColGlobDofs,vector<DofIdx>& DirihletBC);
  virtual PetscErrorCode SetDirihletBC_to_ElementIndicies(
    FieldInterface::FEMethod *fe_method_ptr,
      vector<vector<DofIdx> > &RowGlobDofs,
      vector<vector<DofIdx> > &ColGlobDofs,
      vector<DofIdx>& DirihletBC);

  virtual PetscErrorCode SetDirihletBC_to_ElementIndiciesRow(
    FieldInterface::FEMethod *fe_method_ptr,
    vector<DofIdx>& RowGlobDofs,vector<DofIdx>& DirihletBC);
  virtual PetscErrorCode SetDirihletBC_to_ElementIndiciesCol(
    FieldInterface::FEMethod *fe_method_ptr,
    vector<DofIdx>& ColGlobDofs,vector<DofIdx>& DirihletBC);

  virtual PetscErrorCode SetDirihletBC_to_ElementIndiciesFace(
    FieldInterface::FEMethod *fe_method_ptr,
    vector<DofIdx>& DirihletBC,vector<DofIdx>& FaceNodeGlobalDofs,vector<vector<DofIdx> > &FaceEdgeGlobalDofs,vector<DofIdx> &FaceGlobalDofs);

  virtual PetscErrorCode SetDirihletBC_to_FieldData(FieldInterface::FEMethod *fe_method_ptr,Vec D);
  virtual PetscErrorCode SetDirihletBC_to_MatrixDiagonal(FieldInterface::FEMethod *fe_method_ptr,Mat Aij);
  virtual PetscErrorCode SetDirihletBC_to_RHS(FieldInterface::FEMethod *fe_method_ptr,Vec F);

};

struct CubitDisplacementDirihletBC: public BaseDirihletBC {
  FieldInterface& mField;
  string problem_name;  
  string field_name;

  CubitDisplacementDirihletBC(FieldInterface& _mField,const string _problem_name,const string _field_name); 

  PetscErrorCode ierr;
  ErrorCode rval;

  map<int,Range> bc_map[3];
  map<int,double> bc_map_val[3];

  virtual PetscErrorCode Init();

  PetscErrorCode SetDirihletBC_to_ElementIndiciesRow(
    FieldInterface::FEMethod *fe_method_ptr,vector<vector<DofIdx> > &RowGlobDofs,vector<DofIdx>& DirihletBC);
  PetscErrorCode SetDirihletBC_to_ElementIndiciesCol(
    FieldInterface::FEMethod *fe_method_ptr,vector<vector<DofIdx> > &ColGlobDofs,vector<DofIdx>& DirihletBC);
  PetscErrorCode SetDirihletBC_to_ElementIndicies(
    FieldInterface::FEMethod *fe_method_ptr,vector<vector<DofIdx> > &RowGlobDofs,vector<vector<DofIdx> > &ColGlobDofs,vector<DofIdx>& DirihletBC);
  PetscErrorCode SetDirihletBC_to_ElementIndiciesFace(
    FieldInterface::FEMethod *fe_method_ptr,
    vector<DofIdx>& DirihletBC,vector<DofIdx> &FaceNodeIndices, vector<vector<DofIdx> > &FaceEdgeIndices, vector<DofIdx> &FaceIndices);
  PetscErrorCode SetDirihletBC_to_MatrixDiagonal(FieldInterface::FEMethod *fe_method_ptr,Mat Aij);
  PetscErrorCode SetDirihletBC_to_RHS(FieldInterface::FEMethod *fe_method_ptr,Vec F);
  PetscErrorCode SetDirihletBC_to_FieldData(FieldInterface::FEMethod *fe_method_ptr,Vec D);

};

struct CubitTemperatureDirihletBC: public CubitDisplacementDirihletBC {

  CubitTemperatureDirihletBC(FieldInterface& _mField,const string _problem_name,const string _field_name); 
  PetscErrorCode Init();

};
    
struct CubitDisplacementDirihletBC_ZerosRowsColumns: public CubitDisplacementDirihletBC {
    CubitDisplacementDirihletBC_ZerosRowsColumns(FieldInterface& _mField,const string _problem_name,const string _field_name);
    
    PetscErrorCode SetDirihletBC_to_ElementIndiciesRow(
        FieldInterface::FEMethod *fe_method_ptr,vector<vector<DofIdx> > &RowGlobDofs,vector<DofIdx>& DirihletBC){
        PetscFunctionBegin;
        PetscFunctionReturn(0);
    }
    
    PetscErrorCode SetDirihletBC_to_ElementIndiciesCol(
        FieldInterface::FEMethod *fe_method_ptr,vector<vector<DofIdx> > &ColGlobDofs,vector<DofIdx>& DirihletBC){
        PetscFunctionBegin;
        PetscFunctionReturn(0);
    }
    
    PetscErrorCode SetDirihletBC_to_ElementIndicies(
        FieldInterface::FEMethod *fe_method_ptr,vector<vector<DofIdx> > &RowGlobDofs,vector<vector<DofIdx> > &ColGlobDofs,vector<DofIdx>& DirihletBC){
        PetscFunctionBegin;
        PetscFunctionReturn(0);
    }
    
    PetscErrorCode SetDirihletBC_to_ElementIndiciesFace(
        FieldInterface::FEMethod *fe_method_ptr,
        vector<DofIdx>& DirihletBC,vector<DofIdx> &FaceNodeIndices, vector<vector<DofIdx> > &FaceEdgeIndices, vector<DofIdx> &FaceIndices){
        PetscFunctionBegin;
        PetscFunctionReturn(0);
    }
    
    PetscErrorCode SetDirihletBC_to_MatrixDiagonal(FieldInterface::FEMethod *fe_method_ptr,Mat Aij);
    
};
    

    
}
#endif //__MOABFEMETHOD_DIRIHLETBC_HPP__
