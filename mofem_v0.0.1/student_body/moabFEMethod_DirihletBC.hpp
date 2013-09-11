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
    moabField::FEMethod *fe_method_ptr,string field_name,
    vector<vector<DofIdx> > &RowGlobDofs,vector<vector<DofIdx> > &ColGlobDofs,vector<DofIdx>& DirihletBC) {
    PetscFunctionBegin;
    SETERRQ(PETSC_COMM_SELF,1,"sorry.. you need to tell me what to do");
    NOT_USED(fe_method_ptr);
    NOT_USED(field_name);
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

  virtual PetscErrorCode SetDirihletBC_to_MatrixDiagonal(
    moabField::FEMethod *fe_method_ptr,Mat Aij) {
    PetscFunctionBegin;
    SETERRQ(PETSC_COMM_SELF,1,"sorry.. you need to tell me what to do");
    NOT_USED(fe_method_ptr);
    NOT_USED(Aij);
    PetscFunctionReturn(0);
  }

  virtual PetscErrorCode SetDirihletBC_to_RHS(moabField::FEMethod *fe_method_ptr) {
    PetscFunctionBegin;
    SETERRQ(PETSC_COMM_SELF,1,"sorry.. you need to tell me what to do");
    NOT_USED(fe_method_ptr);
    PetscFunctionReturn(0);
  }

};

struct CubitDirihletBC {
  moabField& mField;
  string field_name;  

  Range nodes[3],edges[3],faces[3];
  Range bc_ents[3];
  CubitDirihletBC(moabField& _mField,string& _field_name): mField(_mField),field_name(_field_name) {};

  PetscErrorCode ierr;
  ErrorCode rval;

  PetscErrorCode Init() {
    PetscFunctionBegin;
    for(_IT_CUBITMESHSETS_FOR_LOOP_(mField,NodeSet|DisplacementSet,it)) {
      displacement_cubit_bc_data mydata;
      ierr = it->get_cubit_bc_data_structure(mydata); CHKERRQ(ierr);
      Range _nodes;
      ierr = it->get_Cubit_msId_entities_by_dimension(mField.get_moab(),0,_nodes,true); CHKERRQ(ierr);
      if(mydata.data.flag1 == 1) nodes[0].insert(_nodes.begin(),_nodes.end());
      if(mydata.data.flag2 == 1) nodes[1].insert(_nodes.begin(),_nodes.end());
      if(mydata.data.flag3 == 1) nodes[2].insert(_nodes.begin(),_nodes.end());
      Range _edges;
      ierr = it->get_Cubit_msId_entities_by_dimension(mField.get_moab(),1,_edges,true); CHKERRQ(ierr);
      rval = mField.get_moab().get_connectivity(_edges,_nodes,true); CHKERR_PETSC(rval);
      if(mydata.data.flag1 == 1) {
	nodes[0].insert(_nodes.begin(),_nodes.end());
	edges[0].insert(_edges.begin(),_edges.end());
      }
      if(mydata.data.flag2 == 1)  {
	nodes[1].insert(_nodes.begin(),_nodes.end());
	edges[1].insert(_edges.begin(),_edges.end());
      }
      if(mydata.data.flag3 == 1)  {
	nodes[2].insert(_nodes.begin(),_nodes.end());
	edges[2].insert(_edges.begin(),_edges.end());
      }
      Range _faces;
      ierr = it->get_Cubit_msId_entities_by_dimension(mField.get_moab(),1,_faces,true); CHKERRQ(ierr);
      rval = mField.get_moab().get_connectivity(_faces,_nodes,true); CHKERR_PETSC(rval);
      ierr = mField.get_moab().get_adjacencies(_faces,1,false,_edges,Interface::UNION); CHKERRQ(ierr);
      if(mydata.data.flag1 == 1) {
	nodes[0].insert(_nodes.begin(),_nodes.end());
	edges[0].insert(_edges.begin(),_edges.end());
	faces[0].insert(_faces.begin(),_faces.end());
      }
      if(mydata.data.flag2 == 1) {
	nodes[1].insert(_nodes.begin(),_nodes.end());
	edges[1].insert(_edges.begin(),_edges.end());
	faces[1].insert(_faces.begin(),_faces.end());
      }
      if(mydata.data.flag3 == 1) {
	nodes[2].insert(_nodes.begin(),_nodes.end());
	edges[2].insert(_edges.begin(),_edges.end());
	faces[2].insert(_faces.begin(),_faces.end());
      }
    }
    for(int ss = 0;ss<3;ss++) {
      bc_ents[ss].insert(nodes[ss].begin(),nodes[ss].end());
      bc_ents[ss].insert(edges[ss].begin(),edges[ss].end());
      bc_ents[ss].insert(faces[ss].begin(),faces[ss].end());
    }
    PetscFunctionReturn(0);
  }

  PetscErrorCode SetDirihletBC_to_ElementIndicies(
    moabField::FEMethod *fe_method_ptr,string field_name,
    vector<vector<DofIdx> > &RowGlobDofs,vector<vector<DofIdx> > &ColGlobDofs,vector<DofIdx>& DirihletBC) {
    PetscFunctionBegin;
    DirihletBC.resize(0);
    for(int ss = 0;ss<3;ss++) {
      Range::iterator siit1 = bc_ents[ss].begin();
      for(;siit1!=bc_ents[ss].end();siit1++) {
	FENumeredDofMoFEMEntity_multiIndex::index<MoABEnt_mi_tag>::type::iterator riit = fe_method_ptr->row_multiIndex->get<MoABEnt_mi_tag>().lower_bound(*siit1);
	FENumeredDofMoFEMEntity_multiIndex::index<MoABEnt_mi_tag>::type::iterator hi_riit = fe_method_ptr->row_multiIndex->get<MoABEnt_mi_tag>().upper_bound(*siit1);
	for(;riit!=hi_riit;riit++) {
	  if(riit->get_name()!=field_name) continue;
	  if(riit->get_dof_rank()!=ss) continue;
	  // if some ranks are selected then we could apply BC in particular direction
	  DirihletBC.push_back(riit->get_petsc_gloabl_dof_idx());
	  for(unsigned int cc = 0;cc<ColGlobDofs.size();cc++) {
	    vector<DofIdx>::iterator it = find(ColGlobDofs[cc].begin(),ColGlobDofs[cc].end(),riit->get_petsc_gloabl_dof_idx());
	    if( it!=ColGlobDofs[cc].end() ) *it = -1; // of idx is set -1 column is not assembled
	  }
	  for(unsigned int rr = 0;rr<RowGlobDofs.size();rr++) {
	    vector<DofIdx>::iterator it = find(RowGlobDofs[rr].begin(),RowGlobDofs[rr].end(),riit->get_petsc_gloabl_dof_idx());
	    if( it!=RowGlobDofs[rr].end() ) *it = -1; // of idx is set -1 row is not assembled
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

};

}

#endif //__MOABFEMETHOD_DIRIHLETBC_HPP__
