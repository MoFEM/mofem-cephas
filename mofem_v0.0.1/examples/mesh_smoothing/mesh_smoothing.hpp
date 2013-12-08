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

#include "FieldInterface.hpp"
#include "FieldCore.hpp"
#include <petscksp.h>

#include "SnesCtx.hpp"
#include "PostProcVertexMethod.hpp"
#include "PostProcDisplacementAndStrainOnRefindedMesh.hpp"

#include "FEMethod_DriverComplexForLazy.hpp"

using namespace MoFEM;

ErrorCode rval;
PetscErrorCode ierr;

static char help[] = "...\n\n";

struct MyMeshSmoothing_ElasticFEMethod_LagnageMultiplaiers: public FEMethod_DriverComplexForLazy_MeshSmoothing  {

  MyMeshSmoothing_ElasticFEMethod_LagnageMultiplaiers(FieldInterface& _mField,BaseDirihletBC *_dirihlet_bc_method_ptr,int _verbose = 0):
    FEMethod_ComplexForLazy_Data(_mField,_dirihlet_bc_method_ptr,_verbose), 
    FEMethod_DriverComplexForLazy_MeshSmoothing(_mField,_dirihlet_bc_method_ptr,0,_verbose) {
    set_qual_ver(1);
  }

};



struct MyMeshSmoothing_ElasticFEMethod: public FEMethod_DriverComplexForLazy_MeshSmoothingProjected {

  MyMeshSmoothing_ElasticFEMethod(FieldInterface& _mField,matPROJ_ctx &_proj_all_ctx,BaseDirihletBC *_dirihlet_bc_method_ptr,int _verbose = 0):
    FEMethod_ComplexForLazy_Data(_mField,_dirihlet_bc_method_ptr,_verbose), 
    FEMethod_DriverComplexForLazy_MeshSmoothingProjected(_mField,_proj_all_ctx,_dirihlet_bc_method_ptr,0,_verbose) {
    set_qual_ver(1);
  }

};

struct materialDirihletBC: public BaseDirihletBC {

  Interface& moab;
  Range &CornersNodes;
  string field_name;
  materialDirihletBC(Interface &_moab,Range& _CornerNodes): moab(_moab),CornersNodes(_CornerNodes),field_name("MESH_NODE_POSITIONS") {}

  PetscErrorCode SetDirihletBC_to_ElementIndicies(
      FieldInterface::FEMethod *fe_method_ptr,vector<vector<DofIdx> > &RowGlob,vector<vector<DofIdx> > &ColGlob,vector<DofIdx>& DirihletBC) {
      PetscFunctionBegin;
      //Dirihlet form SideSet1
      DirihletBC.resize(0);
      Range::iterator siit1 = CornersNodes.begin();
      for(;siit1!=CornersNodes.end();siit1++) {
	  FENumeredDofMoFEMEntity_multiIndex::index<MoABEnt_mi_tag>::type::iterator riit = fe_method_ptr->row_multiIndex->get<MoABEnt_mi_tag>().lower_bound(*siit1);
	  FENumeredDofMoFEMEntity_multiIndex::index<MoABEnt_mi_tag>::type::iterator hi_riit = fe_method_ptr->row_multiIndex->get<MoABEnt_mi_tag>().upper_bound(*siit1);
	  for(;riit!=hi_riit;riit++) {
	    if(riit->get_name()!=field_name) continue;
	    // all fixed
	    DirihletBC.push_back(riit->get_petsc_gloabl_dof_idx());
	    for(unsigned int rr = 0;rr<RowGlob.size();rr++) {
	      vector<DofIdx>::iterator it = find(RowGlob[rr].begin(),RowGlob[rr].end(),riit->get_petsc_gloabl_dof_idx());
	      if( it!=RowGlob[rr].end() ) *it = -1; // of idx is set -1 row is not assembled
	    }
	  }
	  FENumeredDofMoFEMEntity_multiIndex::index<MoABEnt_mi_tag>::type::iterator ciit = fe_method_ptr->col_multiIndex->get<MoABEnt_mi_tag>().lower_bound(*siit1);
	  FENumeredDofMoFEMEntity_multiIndex::index<MoABEnt_mi_tag>::type::iterator hi_ciit = fe_method_ptr->col_multiIndex->get<MoABEnt_mi_tag>().upper_bound(*siit1);
	  for(;ciit!=hi_ciit;ciit++) {
	    if(ciit->get_name()!=field_name) continue;
	    for(unsigned int cc = 0;cc<ColGlob.size();cc++) {
	      vector<DofIdx>::iterator it = find(ColGlob[cc].begin(),ColGlob[cc].end(),ciit->get_petsc_gloabl_dof_idx());
	      if( it!=ColGlob[cc].end() ) *it = -1; // of idx is set -1 column is not assembled
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
    FieldInterface::FEMethod *fe_method_ptr,Mat Aij) {
    PetscFunctionBegin;

    for(_IT_NUMEREDDOFMOFEMENTITY_ROW_BY_LOCIDX_FOR_LOOP_(fe_method_ptr->problem_ptr,dit)) {
      if(find(CornersNodes.begin(),CornersNodes.end(),dit->get_ent()) == CornersNodes.end()) continue;
      ierr = MatSetValue(Aij,dit->get_petsc_gloabl_dof_idx(),dit->get_petsc_gloabl_dof_idx(),1.,INSERT_VALUES); CHKERRQ(ierr);
    }

    PetscFunctionReturn(0);
  }

  PetscErrorCode SetDirihletBC_to_RHS(FieldInterface::FEMethod *fe_method_ptr,Vec F) {
    PetscFunctionBegin;
    for(_IT_NUMEREDDOFMOFEMENTITY_COL_BY_LOCIDX_FOR_LOOP_(fe_method_ptr->problem_ptr,dit)) {
      if(find(CornersNodes.begin(),CornersNodes.end(),dit->get_ent()) == CornersNodes.end()) continue;
      ierr = VecSetValue(F,dit->get_petsc_gloabl_dof_idx(),0.,INSERT_VALUES); CHKERRQ(ierr);
    }
    PetscFunctionReturn(0);
  }

};






