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

#include "Core_dataStructures.hpp"
#include "moabField.hpp"
#include "moabFEMethod_LowLevelStudent.hpp"

#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/matrix_proxy.hpp>
#include <boost/numeric/ublas/io.hpp>

using namespace boost::numeric;

namespace MoFEM {

/** 
 * \brief The student user interface for Dirihlet boundary conditions
 * 
*/
struct BaseDirihletBC {

  BaseDirihletBC() {};

  virtual PetscErrorCode SetDirihletBC_to_ElementIndicies(
    moabField::FEMethod *fe_method_ptr,vector<vector<DofIdx> > &RowGlobDofs,vector<vector<DofIdx> > &ColGlobDofs,vector<DofIdx>& DirihletBC) {
    PetscFunctionBegin;
    SETERRQ(PETSC_COMM_SELF,1,"sorry.. you need to tell me what to do");
    NOT_USED(fe_method_ptr);
    NOT_USED(RowGlobDofs);
    NOT_USED(ColGlobDofs);
    NOT_USED(DirihletBC);
    PetscFunctionReturn(0);
  }

  virtual PetscErrorCode SetDirihletBC_to_ElementIndiciesFace(
    vector<DofIdx>& DirihletBC,vector<DofIdx>& FaceNodeGlobalDofs,vector<vector<DofIdx> > &FaceEdgeGlobalDofs,vector<DofIdx> &FaceGlobalDofs) {
    PetscFunctionBegin;
    SETERRQ(PETSC_COMM_SELF,1,"sorry.. you need to tell me what to do");
    NOT_USED(FaceNodeGlobalDofs);
    NOT_USED(FaceEdgeGlobalDofs);
    NOT_USED(FaceGlobalDofs);
    NOT_USED(DirihletBC);
    PetscFunctionReturn(0);
  }


  virtual PetscErrorCode SetDirihletBC_to_FieldData(moabField::FEMethod *fe_method_ptr,Vec D) {
    PetscFunctionBegin;
    SETERRQ(PETSC_COMM_SELF,1,"sorry.. you need to tell me what to do");
    NOT_USED(fe_method_ptr);
    PetscFunctionReturn(0);
  }

  virtual PetscErrorCode SetDirihletBC_to_MatrixDiagonal(moabField::FEMethod *fe_method_ptr,Mat Aij) {
    PetscFunctionBegin;
    SETERRQ(PETSC_COMM_SELF,1,"sorry.. you need to tell me what to do");
    NOT_USED(fe_method_ptr);
    NOT_USED(Aij);
    PetscFunctionReturn(0);
  }

  virtual PetscErrorCode SetDirihletBC_to_RHS(moabField::FEMethod *fe_method_ptr,Vec F) {
    PetscFunctionBegin;
    SETERRQ(PETSC_COMM_SELF,1,"sorry.. you need to tell me what to do");
    NOT_USED(fe_method_ptr);
    NOT_USED(F);
    PetscFunctionReturn(0);
  }

};

struct CubitDisplacementDirihletBC: public BaseDirihletBC {
  moabField& mField;
  const string problem_name;  
  const string field_name;

  CubitDisplacementDirihletBC(moabField& _mField,const string _problem_name,const string _field_name): 
    mField(_mField),problem_name(_problem_name),field_name(_field_name) {};

  PetscErrorCode ierr;
  ErrorCode rval;

  map<int,Range> bc_map[3];
  map<int,double> bc_map_val[3];

  PetscErrorCode Init() {
    PetscFunctionBegin;
    for(_IT_CUBITMESHSETS_BY_BCDATA_TYPE_FOR_LOOP_(mField,NodeSet|DisplacementSet,it)) {
      ostringstream ss;
      ss << *it << endl;
      displacement_cubit_bc_data mydata;
      ierr = it->get_cubit_bc_data_structure(mydata); CHKERRQ(ierr);
      ss << mydata;
      for(int dim = 0;dim<3;dim++) {
	Range _ents;
	ierr = it->get_Cubit_msId_entities_by_dimension(mField.get_moab(),dim,_ents,true); CHKERRQ(ierr);
	ss << "dim  = " << dim << " nb. ents " << _ents.size() << endl;
	if(dim>1) {
	  Range _edges;
	  ierr = mField.get_moab().get_adjacencies(_ents,1,false,_edges,Interface::UNION); CHKERRQ(ierr);
	  _ents.insert(_edges.begin(),_edges.end());
	  ss << "dim  = " << dim << " nb. edges " << _edges.size() << endl;
	}
	if(dim>0) {
	  Range _nodes;
	  rval = mField.get_moab().get_connectivity(_ents,_nodes,true); CHKERR_PETSC(rval);
	  _ents.insert(_nodes.begin(),_nodes.end());
	  ss << "dim  = " << dim << " nb. nodes " << _nodes.size() << endl;
	}
	if(dim>2) SETERRQ(PETSC_COMM_SELF,1,"not yet implemented");
	if(mydata.data.flag1 == 1) {
	  (bc_map[0])[it->get_msId()].insert(_ents.begin(),_ents.end());
	  (bc_map_val[0])[it->get_msId()] = mydata.data.value1;
	}
	if(mydata.data.flag2 == 1) {
	  (bc_map[1])[it->get_msId()].insert(_ents.begin(),_ents.end());
	  (bc_map_val[1])[it->get_msId()] = mydata.data.value2;
	}
	if(mydata.data.flag3 == 1) {
	  (bc_map[2])[it->get_msId()].insert(_ents.begin(),_ents.end());
	  (bc_map_val[2])[it->get_msId()] = mydata.data.value3;
	}
      }
      ss << endl;
      PetscPrintf(PETSC_COMM_WORLD,ss.str().c_str());
    }
    PetscFunctionReturn(0);
  }

  PetscErrorCode SetDirihletBC_to_ElementIndicies(
    moabField::FEMethod *fe_method_ptr,vector<vector<DofIdx> > &RowGlobDofs,vector<vector<DofIdx> > &ColGlobDofs,vector<DofIdx>& DirihletBC) {
    PetscFunctionBegin;
    DirihletBC.resize(0);
    for(int ss = 0;ss<3;ss++) {
      map<int,Range>::iterator bit = bc_map[ss].begin();
      for(;bit!=bc_map[ss].end();bit++) {
	Range::iterator eit = bit->second.begin();
	for(;eit!=bit->second.end();eit++) {
	  for(_IT_GET_FEROW_DOFS_BY_NAME_AND_ENT_FOR_LOOP_(fe_method_ptr,field_name,*eit,dit)) {
	    if(dit->get_dof_rank()!=ss) continue;
	    // if some ranks are selected then we could apply BC in particular direction
	    DirihletBC.push_back(dit->get_petsc_gloabl_dof_idx());
	    for(unsigned int cc = 0;cc<ColGlobDofs.size();cc++) {
	      vector<DofIdx>::iterator it = find(ColGlobDofs[cc].begin(),ColGlobDofs[cc].end(),dit->get_petsc_gloabl_dof_idx());
	      if( it!=ColGlobDofs[cc].end() ) *it = -1; // of idx is set -1 column is not assembled
	    }
	    for(unsigned int rr = 0;rr<RowGlobDofs.size();rr++) {
	      vector<DofIdx>::iterator it = find(RowGlobDofs[rr].begin(),RowGlobDofs[rr].end(),dit->get_petsc_gloabl_dof_idx());
	      if( it!=RowGlobDofs[rr].end() ) *it = -1; // of idx is set -1 row is not assembled
	    }
	  }
	}
      }
    }
    PetscFunctionReturn(0);
  }

  PetscErrorCode SetDirihletBC_to_ElementIndiciesFace(vector<DofIdx>& DirihletBC,
    vector<DofIdx> &FaceNodeIndices,
    vector<vector<DofIdx> > &FaceEdgeIndices,
    vector<DofIdx> &FaceIndices) {
    PetscFunctionBegin;
    vector<DofIdx>::iterator dit = DirihletBC.begin();
    for(;dit!=DirihletBC.end();dit++) {
      vector<DofIdx>::iterator it = find(FaceNodeIndices.begin(),FaceNodeIndices.end(),*dit);
      if(it!=FaceNodeIndices.end()) *it = -1; // of idx is set -1 row is not assembled
      for(int ee = 0;ee<3;ee++) {
	it = find(FaceEdgeIndices[ee].begin(),FaceEdgeIndices[ee].end(),*dit);
	if(it!=FaceEdgeIndices[ee].end()) *it = -1; // of idx is set -1 row is not assembled
      }
      it = find(FaceIndices.begin(),FaceIndices.end(),*dit);
      if(it!=FaceIndices.end()) *it = -1; // of idx is set -1 row is not assembled
    }
    PetscFunctionReturn(0);
  }

  PetscErrorCode SetDirihletBC_to_MatrixDiagonal(
    moabField::FEMethod *fe_method_ptr,Mat Aij) {
    PetscFunctionBegin;
    ParallelComm* pcomm = ParallelComm::get_pcomm(&mField.get_moab(),MYPCOMM_INDEX);

    for(_IT_NUMEREDDOFMOFEMENTITY_ROW_BY_LOCIDX_FOR_LOOP_(fe_method_ptr->problem_ptr,dit)) {
      if(dit->get_part()!=pcomm->rank()) continue;
      for(int ss = 0;ss<3;ss++) {
	if(dit->get_dof_rank()==ss) {
	  map<int,Range>::iterator bit = bc_map[ss].begin();
	  for(;bit!=bc_map[ss].end();bit++) {
	    if(find(bit->second.begin(),bit->second.end(),dit->get_ent()) == bit->second.end()) continue;
	    ierr = MatSetValue(Aij,dit->get_petsc_gloabl_dof_idx(),dit->get_petsc_gloabl_dof_idx(),1.,INSERT_VALUES); CHKERRQ(ierr);
	  }
	}
      }
    }

    PetscFunctionReturn(0);
  }

  PetscErrorCode SetDirihletBC_to_RHS(moabField::FEMethod *fe_method_ptr,Vec F) {
    PetscFunctionBegin;
    ParallelComm* pcomm = ParallelComm::get_pcomm(&mField.get_moab(),MYPCOMM_INDEX);

    for(_IT_NUMEREDDOFMOFEMENTITY_ROW_BY_LOCIDX_FOR_LOOP_(fe_method_ptr->problem_ptr,dit)) {
      if(dit->get_part()!=pcomm->rank()) continue;
      for(int ss = 0;ss<3;ss++) {
	if(dit->get_dof_rank()==ss) {
	  map<int,Range>::iterator bit = bc_map[ss].begin();
	  for(;bit!=bc_map[ss].end();bit++) {
	    if(find(bit->second.begin(),bit->second.end(),dit->get_ent()) == bit->second.end()) continue;
	    if(dit->get_ent_type()==MBVERTEX) {
	      ierr = VecSetValue(F,dit->get_petsc_gloabl_dof_idx(),(bc_map_val[ss])[bit->first],INSERT_VALUES); CHKERRQ(ierr);
	    } else {
	      ierr = VecSetValue(F,dit->get_petsc_gloabl_dof_idx(),0,INSERT_VALUES); CHKERRQ(ierr);
	    }
	  }
	}
      }
    }

    PetscFunctionReturn(0);
  }

  PetscErrorCode SetDirihletBC_to_FieldData(moabField::FEMethod *fe_method_ptr,Vec D) {
    PetscFunctionBegin;

    ParallelComm* pcomm = ParallelComm::get_pcomm(&mField.get_moab(),MYPCOMM_INDEX);

    ierr = mField.set_local_VecCreateGhost(problem_name,Col,D,INSERT_VALUES,SCATTER_FORWARD); CHKERRQ(ierr);

    for(_IT_NUMEREDDOFMOFEMENTITY_COL_BY_LOCIDX_FOR_LOOP_(fe_method_ptr->problem_ptr,dit)) {
      if(dit->get_part()!=pcomm->rank()) continue;
      if(dit->get_ent_type()!=MBVERTEX) continue;
      for(int ss = 0;ss<3;ss++) {
	if(dit->get_dof_rank()==ss) {
	  map<int,Range>::iterator bit = bc_map[ss].begin();
	  for(;bit!=bc_map[ss].end();bit++) {
	    if(find(bit->second.begin(),bit->second.end(),dit->get_ent()) == bit->second.end()) continue;
	    ierr = VecSetValue(D,dit->get_petsc_gloabl_dof_idx(),(bc_map_val[ss])[bit->first],INSERT_VALUES); CHKERRQ(ierr);
	  }
	}
      }
    }

    ierr = VecAssemblyBegin(D); CHKERRQ(ierr);
    ierr = VecAssemblyEnd(D); CHKERRQ(ierr);

    ierr = VecGhostUpdateBegin(D,INSERT_VALUES,SCATTER_FORWARD); CHKERRQ(ierr);
    ierr = VecGhostUpdateEnd(D,INSERT_VALUES,SCATTER_FORWARD); CHKERRQ(ierr);

    ierr = mField.set_global_VecCreateGhost(problem_name,Col,D,INSERT_VALUES,SCATTER_REVERSE); CHKERRQ(ierr);

    PetscFunctionReturn(0);
  }



};

}

#endif //__MOABFEMETHOD_DIRIHLETBC_HPP__
