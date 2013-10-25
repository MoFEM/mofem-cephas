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


#ifndef __MOABFEMETHOD_DRIVERCOMPLEXFORLAZY_HPP__
#define __MOABFEMETHOD_DRIVERCOMPLEXFORLAZY_HPP__

#include "FieldInterface.hpp"
#include "FieldCore.hpp"
#include "FEMethod_UpLevelStudent.hpp"
#include "DirihletBC.hpp"
#include "cholesky.hpp"
#include <petscksp.h>

#include "SnesCtx.hpp"
#include "FEMethod_ComplexForLazy.hpp"

#include "petscShellMATs_ConstrainsByMarkAinsworth.hpp"
#include "FEMethod_SurfaceConstrains.hpp"

#include "complex_for_lazy.h"

using namespace boost::numeric;

namespace MoFEM {

struct DirihletBCMethod_DriverComplexForLazy: public CubitDisplacementDirihletBC {

  DirihletBCMethod_DriverComplexForLazy(FieldInterface& _mField,const string _problem_name,const string _field_name): 
    CubitDisplacementDirihletBC(_mField,_problem_name,_field_name) {};

};

struct FEMethod_DriverComplexForLazy_Spatial: public FEMethod_ComplexForLazy {

  double *t_val;
  PetscErrorCode set_t_val(double t_val_) {
      PetscFunctionBegin;
      *t_val = t_val_;
      PetscFunctionReturn(0);
  }

  Tag th_t_val;
  FEMethod_DriverComplexForLazy_Spatial(FieldInterface& _mField,BaseDirihletBC *_dirihlet_bc_method_ptr,double _lambda,double _mu,int _verbose = 0): 
  FEMethod_ComplexForLazy_Data(_mField,_dirihlet_bc_method_ptr,_verbose),
  FEMethod_ComplexForLazy(_mField,_dirihlet_bc_method_ptr,FEMethod_ComplexForLazy::spatail_analysis,_lambda,_mu,_verbose) { 
    double def_t_val = 0;
    const EntityHandle root_meshset = mField.get_moab().get_root_set();
    rval = mField.get_moab().tag_get_handle("_LoadFactor_t_val",1,MB_TYPE_DOUBLE,th_t_val,MB_TAG_CREAT|MB_TAG_EXCL|MB_TAG_MESH,&def_t_val); 
    if(rval == MB_ALREADY_ALLOCATED) {
      rval = mField.get_moab().tag_get_by_ptr(th_t_val,&root_meshset,1,(const void**)&t_val); CHKERR_THROW(rval);
    } else {
      CHKERR_THROW(rval);
      rval = mField.get_moab().tag_set_data(th_t_val,&root_meshset,1,&def_t_val); CHKERR_THROW(rval);
      rval = mField.get_moab().tag_get_by_ptr(th_t_val,&root_meshset,1,(const void**)&t_val); CHKERR_THROW(rval);
    }
  };

  PetscLogDouble t1,t2;
  PetscLogDouble v1,v2;

  PetscErrorCode preProcess() {
    PetscFunctionBegin;
    //PetscSynchronizedPrintf(PETSC_COMM_WORLD,"Start Assembly\n");
    ierr = PetscGetTime(&v1); CHKERRQ(ierr);
    ierr = PetscGetCPUTime(&t1); CHKERRQ(ierr);
    switch(snes_ctx) {
      case ctx_SNESNone:
      case ctx_SNESSetFunction: { 
	ierr = VecAssemblyBegin(snes_f); CHKERRQ(ierr);
	ierr = VecAssemblyEnd(snes_f); CHKERRQ(ierr);
      }
      break;
      case ctx_SNESSetJacobian: {
	ierr = MatAssemblyBegin(*snes_B,MAT_FLUSH_ASSEMBLY); CHKERRQ(ierr);
	ierr = MatAssemblyEnd(*snes_B,MAT_FLUSH_ASSEMBLY); CHKERRQ(ierr);
      }
      break;
      default:
	SETERRQ(PETSC_COMM_SELF,1,"not implemented");
    }
    PetscFunctionReturn(0);
  }

  virtual PetscErrorCode CalculateSpatialFint(Vec f) {
    PetscFunctionBegin;

    switch(snes_ctx) {
      case ctx_SNESNone:
      case ctx_SNESSetFunction: { 
        ierr = GetFint(); CHKERRQ(ierr);
	VecSetOption(f,VEC_IGNORE_NEGATIVE_INDICES,PETSC_TRUE); 
	//cerr << "Fint_h " << Fint_h << endl;
	ierr = VecSetValues(f,RowGlobSpatial[i_nodes].size(),&(RowGlobSpatial[i_nodes][0]),&(Fint_h.data()[0]),ADD_VALUES); CHKERRQ(ierr);
	for(int ee = 0;ee<6;ee++) {
	  if(RowGlobSpatial[1+ee].size()>0) {
	    ierr = VecSetValues(f,RowGlobSpatial[1+ee].size(),&(RowGlobSpatial[1+ee][0]),&(Fint_h_edge_data[ee].data()[0]),ADD_VALUES); CHKERRQ(ierr);
	  }
	}
	for(int ff = 0;ff<4;ff++) {
	  if(RowGlobSpatial[1+6+ff].size()>0) {
	    ierr = VecSetValues(f,RowGlobSpatial[1+6+ff].size(),&(RowGlobSpatial[1+6+ff][0]),&(Fint_h_face_data[ff].data()[0]),ADD_VALUES); CHKERRQ(ierr);
	  }
	}
	if(RowGlobSpatial[i_volume].size()>0) {
	  ierr = VecSetValues(f,RowGlobSpatial[i_volume].size(),&(RowGlobSpatial[i_volume][0]),&(Fint_h_volume.data()[0]),ADD_VALUES); CHKERRQ(ierr);
	}
      }
      break;
      default:
	SETERRQ(PETSC_COMM_SELF,1,"not implemented");
    }

    PetscFunctionReturn(0);
  }

  virtual PetscErrorCode CalculateSpatialTangent(Mat B) {
    PetscFunctionBegin;

    switch(snes_ctx) {
      case ctx_SNESNone:
      case ctx_SNESSetJacobian:
	ierr = GetTangent(); CHKERRQ(ierr);
	//cerr << "Khh " << Khh << endl;
	ierr = MatSetValues(B,
	  RowGlobSpatial[i_nodes].size(),&*(RowGlobSpatial[i_nodes].begin()),
	  ColGlobSpatial[i_nodes].size(),&*(ColGlobSpatial[i_nodes].begin()),
	  &*(Khh.data().begin()),ADD_VALUES); CHKERRQ(ierr);
	//
	for(int ee = 0;ee<6;ee++) {
	  ierr = MatSetValues(B,
	    RowGlobSpatial[1+ee].size(),&*(RowGlobSpatial[1+ee].begin()),
	    ColGlobSpatial[i_nodes].size(),&*(ColGlobSpatial[i_nodes].begin()),
	    &*(Kedgeh_data[ee].data().begin()),ADD_VALUES); CHKERRQ(ierr);
	  ierr = MatSetValues(B,
	    RowGlobSpatial[i_nodes].size(),&*(RowGlobSpatial[i_nodes].begin()),
	    ColGlobSpatial[1+ee].size(),&*(ColGlobSpatial[1+ee].begin()),
	    &*(Khedge_data[ee].data().begin()),ADD_VALUES); CHKERRQ(ierr);
	  for(int eee = 0;eee<6;eee++) {
	    ierr = MatSetValues(B,
	      RowGlobSpatial[1+ee].size(),&*(RowGlobSpatial[1+ee].begin()),
	      ColGlobSpatial[1+eee].size(),&*(ColGlobSpatial[1+eee].begin()),
	      &*(Khh_edgeedge_data(ee,eee).data().begin()),ADD_VALUES); CHKERRQ(ierr);
	  }
	  for(int fff = 0;fff<4;fff++) {
	    ierr = MatSetValues(B,
	      RowGlobSpatial[1+ee].size(),&*(RowGlobSpatial[1+ee].begin()),
	      ColGlobSpatial[1+6+fff].size(),&*(ColGlobSpatial[1+6+fff].begin()),
	      &*(Khh_edgeface_data(ee,fff).data().begin()),ADD_VALUES); CHKERRQ(ierr);
	  }
	  ierr = MatSetValues(B,
	    RowGlobSpatial[1+ee].size(),&*(RowGlobSpatial[1+ee].begin()),
	    ColGlobSpatial[i_volume].size(),&*(ColGlobSpatial[i_volume].begin()),
	    &*(Khh_edgevolume_data[ee].data().begin()),ADD_VALUES); CHKERRQ(ierr);
	  ierr = MatSetValues(B,
	    RowGlobSpatial[i_volume].size(),&*(RowGlobSpatial[i_volume].begin()),
	    ColGlobSpatial[1+ee].size(),&*(ColGlobSpatial[1+ee].begin()),
	    &*(Khh_volumeedge_data[ee].data().begin()),ADD_VALUES); CHKERRQ(ierr);
	}
	for(int ff = 0;ff<4;ff++) {
	  ierr = MatSetValues(B,
	    RowGlobSpatial[1+6+ff].size(),&*(RowGlobSpatial[1+6+ff].begin()),
	    ColGlobSpatial[i_nodes].size(),&*(ColGlobSpatial[i_nodes].begin()),
	    &*(Kfaceh_data[ff].data().begin()),ADD_VALUES); CHKERRQ(ierr);
	  ierr = MatSetValues(B,
	    RowGlobSpatial[i_nodes].size(),&*(RowGlobSpatial[i_nodes].begin()),
	    ColGlobSpatial[1+6+ff].size(),&*(ColGlobSpatial[1+6+ff].begin()),
	    &*(Khface_data[ff].data().begin()),ADD_VALUES); CHKERRQ(ierr);
	  for(int eee = 0;eee<6;eee++) {
	    ierr = MatSetValues(B,
	      RowGlobSpatial[1+6+ff].size(),&*(RowGlobSpatial[1+6+ff].begin()),
	      ColGlobSpatial[1+eee].size(),&*(ColGlobSpatial[1+eee].begin()),
	      &*(Khh_faceedge_data(ff,eee).data().begin()),ADD_VALUES); CHKERRQ(ierr);
	  }
	  for(int fff = 0;fff<4;fff++) {
	    ierr = MatSetValues(B,
	      RowGlobSpatial[1+6+ff].size(),&*(RowGlobSpatial[1+6+ff].begin()),
	      ColGlobSpatial[1+6+fff].size(),&*(ColGlobSpatial[1+6+fff].begin()),
	      &*(Khh_faceface_data(ff,fff).data().begin()),ADD_VALUES); CHKERRQ(ierr);
	  }
	  ierr = MatSetValues(B,
	    RowGlobSpatial[1+6+ff].size(),&*(RowGlobSpatial[1+6+ff].begin()),
	    ColGlobSpatial[i_volume].size(),&*(ColGlobSpatial[i_volume].begin()),
	    &*(Khh_facevolume_data[ff].data().begin()),ADD_VALUES); CHKERRQ(ierr);
	  ierr = MatSetValues(B,
	    RowGlobSpatial[i_volume].size(),&*(RowGlobSpatial[i_volume].begin()),
	    ColGlobSpatial[1+6+ff].size(),&*(ColGlobSpatial[1+6+ff].begin()),
	    &*(Khh_volumeface_data[ff].data().begin()),ADD_VALUES); CHKERRQ(ierr);
	}
	ierr = MatSetValues(B,
	  RowGlobSpatial[i_volume].size(),&*(RowGlobSpatial[i_volume].begin()),
	  ColGlobSpatial[i_nodes].size(),&*(ColGlobSpatial[i_nodes].begin()),
	  &*(Kvolumeh.data().begin()),ADD_VALUES); CHKERRQ(ierr);
	ierr = MatSetValues(B,
	  RowGlobSpatial[i_nodes].size(),&*(RowGlobSpatial[i_nodes].begin()),
	  ColGlobSpatial[i_volume].size(),&*(ColGlobSpatial[i_volume].begin()),
	  &*(Khvolume.data().begin()),ADD_VALUES); CHKERRQ(ierr);
	ierr = MatSetValues(B,
	  RowGlobSpatial[i_volume].size(),&*(RowGlobSpatial[i_volume].begin()),
	  ColGlobSpatial[i_volume].size(),&*(ColGlobSpatial[i_volume].begin()),
	  &*(Khh_volumevolume.data().begin()),ADD_VALUES); CHKERRQ(ierr);
	break;
      default:
	SETERRQ(PETSC_COMM_SELF,1,"not implemented");
    }

    PetscFunctionReturn(0);

  }

  virtual PetscErrorCode CaluclateSpatialFext(Vec f,double *t,Range& NeumannSideSet) {
    PetscFunctionBegin;

    SideNumber_multiIndex& side_table = const_cast<SideNumber_multiIndex&>(fe_ent_ptr->get_side_number_table());
    SideNumber_multiIndex::nth_index<1>::type::iterator siit = side_table.get<1>().lower_bound(boost::make_tuple(MBTRI,0));
    SideNumber_multiIndex::nth_index<1>::type::iterator hi_siit = side_table.get<1>().upper_bound(boost::make_tuple(MBTRI,4));

    switch(snes_ctx) {
      case ctx_SNESNone:
      case ctx_SNESSetFunction: { 
	for(;siit!=hi_siit;siit++) {
	  VecSetOption(f,VEC_IGNORE_NEGATIVE_INDICES,PETSC_TRUE); 
	  Range::iterator fit = find(NeumannSideSet.begin(),NeumannSideSet.end(),siit->ent);
	  if(fit==NeumannSideSet.end()) continue;
	  ierr = GetFaceIndicesAndData(siit->ent); CHKERRQ(ierr);
	  ierr = GetFExt(siit->ent,t,NULL,NULL); CHKERRQ(ierr);
	  //cerr << "FExt " << FExt << endl;
	  //cerr << "FaceNodeIndices.size() " << FaceNodeIndices.size() << endl;
	  ierr = dirihlet_bc_method_ptr->SetDirihletBC_to_ElementIndiciesFace(DirihletBC,FaceNodeIndices,FaceEdgeIndices_data,FaceIndices); CHKERRQ(ierr);
	  ierr = VecSetValues(f,FaceNodeIndices.size(),&(FaceNodeIndices[0]),&*FExt.data().begin(),ADD_VALUES); CHKERRQ(ierr);
	  for(int ee = 0;ee<3;ee++) {
	    if(FaceEdgeIndices_data[ee].size()>0) {
	      ierr = VecSetValues(f,FaceEdgeIndices_data[ee].size(),&(FaceEdgeIndices_data[ee][0]),&*FExt_edge_data[ee].data().begin(),ADD_VALUES); CHKERRQ(ierr);
	    }
	  }
	  if(FaceIndices.size()>0) {
	    ierr = VecSetValues(f,FaceIndices.size(),&(FaceIndices[0]),&*FExt_face.data().begin(),ADD_VALUES); CHKERRQ(ierr);
	  }
	}
      }
      break;
      default:
	SETERRQ(PETSC_COMM_SELF,1,"not implemented");
    }

    PetscFunctionReturn(0);
  }

  virtual PetscErrorCode CalculateSpatialTangentExt(Mat B,double *t,Range& NeumannSideSet) {
    PetscFunctionBegin;
    if(get_PhysicalEquationNumber()==hooke) PetscFunctionReturn(0);
  
    SideNumber_multiIndex& side_table = const_cast<SideNumber_multiIndex&>(fe_ent_ptr->get_side_number_table());
    SideNumber_multiIndex::nth_index<1>::type::iterator siit = side_table.get<1>().lower_bound(boost::make_tuple(MBTRI,0));
    SideNumber_multiIndex::nth_index<1>::type::iterator hi_siit = side_table.get<1>().upper_bound(boost::make_tuple(MBTRI,4));

    switch(snes_ctx) {
      case ctx_SNESSetJacobian:
	for(;siit!=hi_siit;siit++) {
	  Range::iterator fit = find(NeumannSideSet.begin(),NeumannSideSet.end(),siit->ent);
	  if(fit==NeumannSideSet.end()) continue;
	  ierr = GetFaceIndicesAndData(siit->ent); CHKERRQ(ierr);
	  ierr = GetTangentExt(siit->ent,t,NULL,NULL); CHKERRQ(ierr);
	  ierr = dirihlet_bc_method_ptr->SetDirihletBC_to_ElementIndiciesFace(DirihletBC,FaceNodeIndices,FaceEdgeIndices_data,FaceIndices); CHKERRQ(ierr);
	  ierr = MatSetValues(B,
	    FaceNodeIndices.size(),&(FaceNodeIndices[0]),FaceNodeIndices.size(),&(FaceNodeIndices[0]),
	    &*(KExt_hh.data().begin()),ADD_VALUES); CHKERRQ(ierr);
	  for(int ee = 0;ee<3;ee++) {
	    if(FaceNodeIndices.size()==0) continue;
	    ierr = MatSetValues(B,
	      FaceEdgeIndices_data[ee].size(),&(FaceEdgeIndices_data[ee][0]),
	      FaceNodeIndices.size(),&(FaceNodeIndices[0]),
	      &*(KExt_edgeh_data[ee].data().begin()),ADD_VALUES); CHKERRQ(ierr);
	    ierr = MatSetValues(B,
	      FaceNodeIndices.size(),&(FaceNodeIndices[0]),
	      FaceEdgeIndices_data[ee].size(),&(FaceEdgeIndices_data[ee][0]),
	      &*(KExt_hedge_data[ee].data().begin()),ADD_VALUES); CHKERRQ(ierr);
	    for(int eee = 0;eee<3;eee++) {
	      ierr = MatSetValues(B,
		FaceEdgeIndices_data[ee].size(),&(FaceEdgeIndices_data[ee][0]),
		FaceEdgeIndices_data[eee].size(),&(FaceEdgeIndices_data[eee][0]),
		&*(KExt_edgeedge_data(ee,eee).data().begin()),ADD_VALUES); CHKERRQ(ierr);
	    }
	    ierr = MatSetValues(B,
	      FaceIndices.size(),&(FaceIndices[0]),
	      FaceEdgeIndices_data[ee].size(),&(FaceEdgeIndices_data[ee][0]),
	      &*(KExt_faceedge_data[ee].data().begin()),ADD_VALUES); CHKERRQ(ierr);
	    ierr = MatSetValues(B,
	      FaceEdgeIndices_data[ee].size(),&(FaceEdgeIndices_data[ee][0]),
	      FaceIndices.size(),&(FaceIndices[0]),
	      &*(KExt_edgeface_data[ee].data().begin()),ADD_VALUES); CHKERRQ(ierr);
	  }
	  ierr = MatSetValues(B,
	    FaceNodeIndices.size(),&(FaceNodeIndices[0]),FaceIndices.size(),&(FaceIndices[0]),
	    &*(KExt_hface.data().begin()),ADD_VALUES); CHKERRQ(ierr);
	  ierr = MatSetValues(B,
	    FaceIndices.size(),&(FaceIndices[0]),FaceNodeIndices.size(),&(FaceNodeIndices[0]),
	    &*(KExt_faceh.data().begin()),ADD_VALUES); CHKERRQ(ierr);
	  ierr = MatSetValues(B,
	    FaceIndices.size(),&(FaceIndices[0]),FaceIndices.size(),&(FaceIndices[0]),
	    &*(KExt_faceface.data().begin()),ADD_VALUES); CHKERRQ(ierr);
	}
	break;
      default:
	SETERRQ(PETSC_COMM_SELF,1,"not implemented");
    }


    PetscFunctionReturn(0);
  }

  vector<DofIdx> DirihletBC;
  PetscErrorCode operator()() {
    PetscFunctionBegin;

    ierr = OpComplexForLazyStart(); CHKERRQ(ierr);
    ierr = GetIndicesSpatial(); CHKERRQ(ierr);

    ierr = dirihlet_bc_method_ptr->SetDirihletBC_to_ElementIndicies(this,RowGlobSpatial,ColGlobSpatial,DirihletBC); CHKERRQ(ierr);

    switch(snes_ctx) {
      case ctx_SNESSetFunction: { 
	ierr = CalculateSpatialFint(snes_f); CHKERRQ(ierr);
      }
      break;
      case ctx_SNESSetJacobian:
	ierr = CalculateSpatialTangent(*snes_B); CHKERRQ(ierr);
	break;
      default:
	SETERRQ(PETSC_COMM_SELF,1,"not implemented");
    }

    for(_IT_CUBITMESHSETS_BY_BCDATA_TYPE_FOR_LOOP_(mField,SideSet|PressureSet,it)) {

      pressure_cubit_bc_data mydata;
      ierr = it->get_cubit_bc_data_structure(mydata); CHKERRQ(ierr);
      /*ostringstream ss;
      ss << *it << endl;
      ss << mydata;
      PetscPrintf(PETSC_COMM_WORLD,ss.str().c_str());*/

      double t_val_ = *(this->t_val)*mydata.data.value1;
      double t[] = { 0,0,-t_val_, 0,0,-t_val_, 0,0,-t_val_ };

      Range NeumannSideSet;
      ierr = it->get_Cubit_msId_entities_by_dimension(mField.get_moab(),2,NeumannSideSet,true); CHKERRQ(ierr);

      switch(snes_ctx) {
	case ctx_SNESSetFunction: { 
	  ierr = CaluclateSpatialFext(snes_f,t,NeumannSideSet); CHKERRQ(ierr);
	}
	break;
	case ctx_SNESSetJacobian:
	  ierr = CalculateSpatialTangentExt(*snes_B,t,NeumannSideSet); CHKERRQ(ierr);
	  break;
	default:
	  SETERRQ(PETSC_COMM_SELF,1,"not implemented");
      }

    }

    PetscFunctionReturn(0);
  }

  PetscErrorCode postProcess() {
    PetscFunctionBegin;
    switch(snes_ctx) {
      case ctx_SNESSetFunction: { 
	ierr = VecAssemblyBegin(snes_f); CHKERRQ(ierr);
	ierr = VecAssemblyEnd(snes_f); CHKERRQ(ierr);
	ierr = dirihlet_bc_method_ptr->SetDirihletBC_to_RHS(this,snes_f); CHKERRQ(ierr);
      }
      break;
      case ctx_SNESSetJacobian: {
	ierr = MatAssemblyBegin(*snes_B,MAT_FLUSH_ASSEMBLY); CHKERRQ(ierr);
	ierr = MatAssemblyEnd(*snes_B,MAT_FLUSH_ASSEMBLY); CHKERRQ(ierr);
	ierr = dirihlet_bc_method_ptr->SetDirihletBC_to_MatrixDiagonal(this,*snes_B); CHKERRQ(ierr);
      }
      break;
      default:
	SETERRQ(PETSC_COMM_SELF,1,"not implemented");
    }

    ierr = PetscGetTime(&v2); CHKERRQ(ierr);
    ierr = PetscGetCPUTime(&t2); CHKERRQ(ierr);
    //PetscSynchronizedPrintf(PETSC_COMM_WORLD,"End Assembly: Rank %d Time = %f CPU Time = %f\n",pcomm->rank(),v2-v1,t2-t1);
    PetscFunctionReturn(0);
  }

};

struct FEMethod_DriverComplexForLazy_Material: public FEMethod_DriverComplexForLazy_Spatial {

  FEMethod_DriverComplexForLazy_Material(FieldInterface& _mField,BaseDirihletBC *_dirihlet_bc_method_ptr,double _lambda,double _mu,int _verbose = 0):
    FEMethod_ComplexForLazy_Data(_mField,_dirihlet_bc_method_ptr,_verbose), 
    FEMethod_DriverComplexForLazy_Spatial(_mField,_dirihlet_bc_method_ptr,_lambda,_mu,_verbose) {
      type_of_analysis = material_analysis;
    }

  virtual PetscErrorCode CalculateMaterialTangent(Mat B) {
    PetscFunctionBegin;

    switch(snes_ctx) {
      case ctx_SNESNone:
      case ctx_SNESSetJacobian:
	ierr = GetTangent(); CHKERRQ(ierr);
	ierr = MatSetValues(B,
	  RowGlobMaterial[i_nodes].size(),&*(RowGlobMaterial[i_nodes].begin()),
	  ColGlobMaterial[i_nodes].size(),&*(ColGlobMaterial[i_nodes].begin()),
	  &*(KHH.data().begin()),ADD_VALUES); CHKERRQ(ierr);
	break;
      default:
	SETERRQ(PETSC_COMM_SELF,1,"not implemented");
    }

    PetscFunctionReturn(0);
  }

  virtual PetscErrorCode CalculateMaterialFint(Vec f) {
    PetscFunctionBegin;

    switch(snes_ctx) {
      case ctx_SNESNone:
      case ctx_SNESSetFunction: { 
        ierr = GetFint(); CHKERRQ(ierr);
	ierr = VecSetOption(f,VEC_IGNORE_NEGATIVE_INDICES,PETSC_TRUE);  CHKERRQ(ierr);
	//cerr << "Fint_h " << Fint_h << endl;
	ierr = VecSetValues(f,RowGlobMaterial[i_nodes].size(),&(RowGlobMaterial[i_nodes][0]),&(Fint_H.data()[0]),ADD_VALUES); CHKERRQ(ierr);
      }
      break;
      default:
	SETERRQ(PETSC_COMM_SELF,1,"not implemented");
    }

    PetscFunctionReturn(0);
  }

  virtual PetscErrorCode CaluclateMaterialFext(Vec f,double *t,Range& NeumannSideSet) {
    PetscFunctionBegin;

    SideNumber_multiIndex& side_table = const_cast<SideNumber_multiIndex&>(fe_ent_ptr->get_side_number_table());
    SideNumber_multiIndex::nth_index<1>::type::iterator siit = side_table.get<1>().lower_bound(boost::make_tuple(MBTRI,0));
    SideNumber_multiIndex::nth_index<1>::type::iterator hi_siit = side_table.get<1>().upper_bound(boost::make_tuple(MBTRI,4));

    vector<vector<DofIdx> > dummy1;
    vector<DofIdx> dummy2;

    switch(snes_ctx) {
      case ctx_SNESNone:
      case ctx_SNESSetFunction: { 
	for(;siit!=hi_siit;siit++) {
	  VecSetOption(f,VEC_IGNORE_NEGATIVE_INDICES,PETSC_TRUE); 
	  Range::iterator fit = find(NeumannSideSet.begin(),NeumannSideSet.end(),siit->ent);
	  if(fit==NeumannSideSet.end()) continue;
	  ierr = GetFaceIndicesAndData_Material(siit->ent); CHKERRQ(ierr);
	  ierr = GetFExt_Material(siit->ent,t,NULL,NULL); CHKERRQ(ierr);
	  ierr = dirihlet_bc_method_ptr->SetDirihletBC_to_ElementIndiciesFace(DirihletBC,FaceNodeIndices_Material,dummy1,dummy2); CHKERRQ(ierr);
	  ierr = VecSetValues(f,FaceNodeIndices_Material.size(),&(FaceNodeIndices_Material[0]),&*FExt_Material.data().begin(),ADD_VALUES); CHKERRQ(ierr);
	}
      }
      break;
      default:
	SETERRQ(PETSC_COMM_SELF,1,"not implemented");
    }

    PetscFunctionReturn(0);
  }

  virtual PetscErrorCode CalculateMaterialTangentExt(Mat B,double *t,Range& NeumannSideSet) {
    PetscFunctionBegin;
    if(get_PhysicalEquationNumber()==hooke) PetscFunctionReturn(0);
  
    SideNumber_multiIndex& side_table = const_cast<SideNumber_multiIndex&>(fe_ent_ptr->get_side_number_table());
    SideNumber_multiIndex::nth_index<1>::type::iterator siit = side_table.get<1>().lower_bound(boost::make_tuple(MBTRI,0));
    SideNumber_multiIndex::nth_index<1>::type::iterator hi_siit = side_table.get<1>().upper_bound(boost::make_tuple(MBTRI,4));

    vector<vector<DofIdx> > dummy1;
    vector<DofIdx> dummy2;

    switch(snes_ctx) {
      case ctx_SNESSetJacobian:
	for(;siit!=hi_siit;siit++) {
	  Range::iterator fit = find(NeumannSideSet.begin(),NeumannSideSet.end(),siit->ent);
	  if(fit==NeumannSideSet.end()) continue;
	  ierr = GetFaceIndicesAndData_Material(siit->ent); CHKERRQ(ierr);
	  ierr = GetTangentExt_Material(siit->ent,t,NULL,NULL); CHKERRQ(ierr);
	  ierr = dirihlet_bc_method_ptr->SetDirihletBC_to_ElementIndiciesFace(DirihletBC,FaceNodeIndices_Material,dummy1,dummy2); CHKERRQ(ierr);
	  ierr = MatSetValues(B,
	    FaceNodeIndices.size(),&(FaceNodeIndices_Material[0]),FaceNodeIndices_Material.size(),&(FaceNodeIndices_Material[0]),
	    &*(KExt_HH_Material.data().begin()),ADD_VALUES); CHKERRQ(ierr);
	}
	break;
      default:
	SETERRQ(PETSC_COMM_SELF,1,"not implemented");
    }

    PetscFunctionReturn(0);
  }


  PetscErrorCode operator()(Mat *B) {
    PetscFunctionBegin;

    ierr = OpComplexForLazyStart(); CHKERRQ(ierr);
    ierr = GetIndicesMaterial(); CHKERRQ(ierr);
    ierr = GetData(dofs_x_edge_data,dofs_x_edge,
      dofs_x_face_data,dofs_x_face,
      dofs_x_volume,dofs_x,
      spatial_field_name); CHKERRQ(ierr);

    ierr = dirihlet_bc_method_ptr->SetDirihletBC_to_ElementIndicies(this,RowGlobMaterial,ColGlobMaterial,DirihletBC); CHKERRQ(ierr);

    switch(snes_ctx) {
      case ctx_SNESSetFunction: { 
	ierr = CalculateMaterialFint(snes_f); CHKERRQ(ierr);
      }
      break;
      case ctx_SNESSetJacobian:
	ierr = CalculateMaterialTangent(*B); CHKERRQ(ierr);
	break;
      default:
	SETERRQ(PETSC_COMM_SELF,1,"not implemented");
    }

    for(_IT_CUBITMESHSETS_BY_BCDATA_TYPE_FOR_LOOP_(mField,SideSet|PressureSet,it)) {

      pressure_cubit_bc_data mydata;
      ierr = it->get_cubit_bc_data_structure(mydata); CHKERRQ(ierr);
      /*ostringstream ss;
      ss << *it << endl;
      ss << mydata;
      PetscPrintf(PETSC_COMM_WORLD,ss.str().c_str());*/

      double t_val_ = *(this->t_val)*mydata.data.value1;
      double t[] = { 0,0,-t_val_, 0,0,-t_val_, 0,0,-t_val_ };

      Range NeumannSideSet;
      ierr = it->get_Cubit_msId_entities_by_dimension(mField.get_moab(),2,NeumannSideSet,true); CHKERRQ(ierr);

      switch(snes_ctx) {
	case ctx_SNESSetFunction: { 
	  ierr = CaluclateMaterialFext(snes_f,t,NeumannSideSet) ; CHKERRQ(ierr);
	}
	break;
	case ctx_SNESSetJacobian:
	  ierr = CalculateMaterialTangentExt(*B,t,NeumannSideSet); CHKERRQ(ierr);
	  break;
	default:
	  SETERRQ(PETSC_COMM_SELF,1,"not implemented");
      }

    }


    PetscFunctionReturn(0);
  }

  PetscErrorCode operator()() {
    PetscFunctionBegin;
    ierr = operator()(snes_B); CHKERRQ(ierr);
    PetscFunctionReturn(0);
  }

};

struct FEMethod_DriverComplexForLazy_MeshSmoothing: public FEMethod_DriverComplexForLazy_Material {

  FEMethod_DriverComplexForLazy_MeshSmoothing(FieldInterface& _mField,BaseDirihletBC *_dirihlet_bc_method_ptr,int _verbose = 0):
    FEMethod_ComplexForLazy_Data(_mField,_dirihlet_bc_method_ptr,_verbose), 
    FEMethod_DriverComplexForLazy_Material(_mField,_dirihlet_bc_method_ptr,0,0,_verbose) {
      type_of_analysis = mesh_quality_analysis;

      g_NTET.resize(4*1);
      ShapeMBTET(&g_NTET[0],G_TET_X1,G_TET_Y1,G_TET_Z1,1);
      g_TET_W = G_TET_W1;

    }

  virtual PetscErrorCode CalculateMeshSmoothingTangent(Mat B) {
    PetscFunctionBegin;

    switch(snes_ctx) {
      case ctx_SNESNone:
      case ctx_SNESSetJacobian:
	ierr = GetTangent(); CHKERRQ(ierr);
	ierr = MatSetValues(B,
	  RowGlobMaterial[i_nodes].size(),&*(RowGlobMaterial[i_nodes].begin()),
	  ColGlobMaterial[i_nodes].size(),&*(ColGlobMaterial[i_nodes].begin()),
	  &*(KHH.data().begin()),ADD_VALUES); CHKERRQ(ierr);
	break;
      default:
	SETERRQ(PETSC_COMM_SELF,1,"not implemented");
    }

    PetscFunctionReturn(0);
  }

  virtual PetscErrorCode CalculateMeshSmoothingFint(Vec f) {
    PetscFunctionBegin;

    switch(snes_ctx) {
      case ctx_SNESNone:
      case ctx_SNESSetFunction: { 
        ierr = GetFint(); CHKERRQ(ierr);
	ierr = VecSetOption(f,VEC_IGNORE_NEGATIVE_INDICES,PETSC_TRUE);  CHKERRQ(ierr);
	//cerr << "Fint_h " << Fint_h << endl;
	ierr = VecSetValues(f,RowGlobMaterial[i_nodes].size(),&(RowGlobMaterial[i_nodes][0]),&(Fint_H.data()[0]),ADD_VALUES); CHKERRQ(ierr);
      }
      break;
      default:
	SETERRQ(PETSC_COMM_SELF,1,"not implemented");
    }

    PetscFunctionReturn(0);
  }


  PetscErrorCode operator()(Mat *B) {
    PetscFunctionBegin;

    ierr = OpComplexForLazyStart(); CHKERRQ(ierr);
    ierr = GetIndicesMaterial(); CHKERRQ(ierr);

    ierr = dirihlet_bc_method_ptr->SetDirihletBC_to_ElementIndicies(this,RowGlobMaterial,ColGlobMaterial,DirihletBC); CHKERRQ(ierr);
  
    switch(snes_ctx) {
      case ctx_SNESSetFunction:  
	ierr = CalculateMeshSmoothingFint(snes_f); CHKERRQ(ierr);
      break;
      case ctx_SNESSetJacobian:
	ierr = CalculateMeshSmoothingTangent(*B); CHKERRQ(ierr);
	break;
      default:
	SETERRQ(PETSC_COMM_SELF,1,"not implemented");
    }

    PetscFunctionReturn(0);
  }

  PetscErrorCode operator()() {
    PetscFunctionBegin;
    ierr = operator()(snes_B); CHKERRQ(ierr);
    PetscFunctionReturn(0);
  }

};

struct FEMethod_DriverComplexForLazy_Projected: public virtual FEMethod_ComplexForLazy_Data {

  matPROJ_ctx& proj_all_ctx;
  bool init;
  double alpha;

  string problem_name;

  FEMethod_DriverComplexForLazy_Projected(FieldInterface& _mField,matPROJ_ctx &_proj_all_ctx,BaseDirihletBC *_dirihlet_bc_method_ptr,string _problem_name): 
    FEMethod_ComplexForLazy_Data(_mField,_dirihlet_bc_method_ptr), 
    proj_all_ctx(_proj_all_ctx),init(true),problem_name(_problem_name) {};

  Range CornersEdges,CornersNodes,SurfacesFaces;
  C_SURFACE_FEMethod *CFE_SURFACE;
  g_SURFACE_FEMethod *gFE_SURFACE;
  C_CORNER_FEMethod *CFE_CORNER;
  g_CORNER_FEMethod *gFE_CORNER;
  //CRACK
  C_SURFACE_FEMethod *CFE_CRACK_SURFACE;
  g_SURFACE_FEMethod *gFE_CRACK_SURFACE;

  PetscErrorCode set_local_VecCreateGhost_for_ConstrainsProblem(Vec x) {
    PetscFunctionBegin;
    Vec _x_;
    ierr = mField.VecCreateGhost("C_ALL_MATRIX",Col,&_x_); CHKERRQ(ierr);
    ierr = VecScatterBegin(proj_all_ctx.scatter,x,_x_,INSERT_VALUES,SCATTER_FORWARD); CHKERRQ(ierr);
    ierr = VecScatterEnd(proj_all_ctx.scatter,x,_x_,INSERT_VALUES,SCATTER_FORWARD); CHKERRQ(ierr);
    ierr = VecGhostUpdateBegin(_x_,INSERT_VALUES,SCATTER_FORWARD); CHKERRQ(ierr);
    ierr = VecGhostUpdateEnd(_x_,INSERT_VALUES,SCATTER_FORWARD); CHKERRQ(ierr);
    ierr = mField.set_local_VecCreateGhost("C_ALL_MATRIX",Col,_x_,INSERT_VALUES,SCATTER_REVERSE); CHKERRQ(ierr);
    ierr = VecDestroy(&_x_); CHKERRQ(ierr);
    PetscFunctionReturn(0);
  }

  PetscErrorCode _preProcess_ctx_SNESSetFunction(Vec x,Vec f) {
    PetscFunctionBegin;

    PetscBool flg;
    ierr = PetscOptionsGetReal("","-my_penalty",&alpha,&flg); CHKERRQ(ierr);
    if(flg != PETSC_TRUE) {
      SETERRQ(PETSC_COMM_SELF,1,"*** ERROR -my_penalty need to be given");
    }

    if(init) {
	init = false;
	ierr = mField.get_Cubit_msId_entities_by_dimension(100,SideSet,1,CornersEdges,true); CHKERRQ(ierr);
	ierr = mField.get_Cubit_msId_entities_by_dimension(101,NodeSet,0,CornersNodes,true); CHKERRQ(ierr);
	ierr = mField.get_Cubit_msId_entities_by_dimension(102,SideSet,2,SurfacesFaces,true); CHKERRQ(ierr);
	ierr = mField.set_global_VecCreateGhost(problem_name,Col,x,INSERT_VALUES,SCATTER_REVERSE); CHKERRQ(ierr);

	Range CornersEdgesNodes;
	rval = moab.get_connectivity(CornersEdges,CornersEdgesNodes,true); CHKERR_PETSC(rval);
	CornersNodes.insert(CornersEdgesNodes.begin(),CornersEdgesNodes.end());
	CFE_SURFACE = new C_SURFACE_FEMethod(moab,SurfacesFaces,proj_all_ctx.C);
	gFE_SURFACE = new g_SURFACE_FEMethod(moab,SurfacesFaces,proj_all_ctx.g);
	CFE_CORNER = new C_CORNER_FEMethod(moab,CornersNodes,proj_all_ctx.C);
	gFE_CORNER = new g_CORNER_FEMethod(moab,CornersNodes,proj_all_ctx.g);
	//CRACK
	CFE_CRACK_SURFACE = new C_SURFACE_FEMethod(moab,SurfacesFaces,proj_all_ctx.C,"LAMBDA_CRACK_SURFACE",0);
	gFE_CRACK_SURFACE = new g_SURFACE_FEMethod(moab,SurfacesFaces,proj_all_ctx.g,"LAMBDA_CRACK_SURFACE",0);

	ierr = MatSetOption(proj_all_ctx.C,MAT_NEW_NONZERO_LOCATION_ERR,PETSC_TRUE); CHKERRQ(ierr);
	ierr = MatSetOption(proj_all_ctx.C,MAT_NEW_NONZERO_ALLOCATION_ERR,PETSC_TRUE); CHKERRQ(ierr);

	ierr = MatZeroEntries(proj_all_ctx.C); CHKERRQ(ierr);
	ierr = mField.loop_finite_elements("C_ALL_MATRIX","C_SURFACE_ELEM",*CFE_SURFACE);  CHKERRQ(ierr);
	ierr = mField.loop_finite_elements("C_ALL_MATRIX","C_CRACK_SURFACE_ELEM",*CFE_CRACK_SURFACE);  CHKERRQ(ierr);
	ierr = MatAssemblyBegin(proj_all_ctx.C,MAT_FLUSH_ASSEMBLY); CHKERRQ(ierr);
	ierr = MatAssemblyEnd(proj_all_ctx.C,MAT_FLUSH_ASSEMBLY); CHKERRQ(ierr);
	ierr = mField.loop_finite_elements("C_ALL_MATRIX","C_CORNER_ELEM",*CFE_CORNER);  CHKERRQ(ierr);
	ierr = MatAssemblyBegin(proj_all_ctx.C,MAT_FINAL_ASSEMBLY); CHKERRQ(ierr);
	ierr = MatAssemblyEnd(proj_all_ctx.C,MAT_FINAL_ASSEMBLY); CHKERRQ(ierr);

	/*{
	  MatView(proj_all_ctx.C,PETSC_VIEWER_DRAW_WORLD);
	  std::string wait;
	  std::cin >> wait;
	}*/

	ierr = proj_all_ctx.InitQorP(f); CHKERRQ(ierr);
	ierr = proj_all_ctx.InitQTKQ(); CHKERRQ(ierr);

	//ierr = proj_all_ctx.RecalculateCTandCCT(); CHKERRQ(ierr);
	//ierr = proj_all_ctx.RecalulateCTC(); CHKERRQ(ierr);

    }  else {
	ierr = set_local_VecCreateGhost_for_ConstrainsProblem(x); CHKERRQ(ierr);
    }

    ierr = VecZeroEntries(proj_all_ctx.g); CHKERRQ(ierr);
    ierr = VecGhostUpdateBegin(proj_all_ctx.g,INSERT_VALUES,SCATTER_FORWARD); CHKERRQ(ierr);
    ierr = VecGhostUpdateEnd(proj_all_ctx.g,INSERT_VALUES,SCATTER_FORWARD); CHKERRQ(ierr);

    ierr = mField.loop_finite_elements("C_ALL_MATRIX","C_SURFACE_ELEM",*gFE_SURFACE);  CHKERRQ(ierr);
    ierr = mField.loop_finite_elements("C_ALL_MATRIX","C_CRACK_SURFACE_ELEM",*gFE_CRACK_SURFACE);  CHKERRQ(ierr);
    ierr = VecGhostUpdateBegin(proj_all_ctx.g,ADD_VALUES,SCATTER_REVERSE); CHKERRQ(ierr);
    ierr = VecGhostUpdateEnd(proj_all_ctx.g,ADD_VALUES,SCATTER_REVERSE); CHKERRQ(ierr);
    ierr = VecAssemblyBegin(proj_all_ctx.g); CHKERRQ(ierr);
    ierr = VecAssemblyEnd(proj_all_ctx.g); CHKERRQ(ierr);

    ierr = mField.loop_finite_elements("C_ALL_MATRIX","C_CORNER_ELEM",*gFE_CORNER);  CHKERRQ(ierr);
    ierr = VecGhostUpdateBegin(proj_all_ctx.g,INSERT_VALUES,SCATTER_REVERSE); CHKERRQ(ierr);
    ierr = VecGhostUpdateEnd(proj_all_ctx.g,INSERT_VALUES,SCATTER_REVERSE); CHKERRQ(ierr);
    ierr = VecAssemblyBegin(proj_all_ctx.g); CHKERRQ(ierr);
    ierr = VecAssemblyEnd(proj_all_ctx.g); CHKERRQ(ierr);

    PetscReal g_nrm2;
    ierr = VecNorm(proj_all_ctx.g, NORM_2,&g_nrm2); CHKERRQ(ierr);
    PetscPrintf(PETSC_COMM_WORLD,"\t g_nrm2 = %6.4e\n",g_nrm2);

    PetscFunctionReturn(0);
  }

  PetscErrorCode _postProcess_ctx_SNESSetFunction(Vec f) {
    PetscFunctionBegin;

	ierr = VecGhostUpdateBegin(f,ADD_VALUES,SCATTER_REVERSE); CHKERRQ(ierr);
	ierr = VecGhostUpdateEnd(f,ADD_VALUES,SCATTER_REVERSE); CHKERRQ(ierr);
	ierr = VecAssemblyBegin(f); CHKERRQ(ierr);
	ierr = VecAssemblyEnd(f); CHKERRQ(ierr);
	int M,N,m,n;
	ierr = MatGetSize(proj_all_ctx.K,&M,&N); CHKERRQ(ierr);
	ierr = MatGetLocalSize(proj_all_ctx.K,&m,&n); CHKERRQ(ierr);
	Mat Q;
	ierr = MatCreateShell(PETSC_COMM_WORLD,m,n,M,N,&proj_all_ctx,&Q); CHKERRQ(ierr);
	ierr = MatShellSetOperation(Q,MATOP_MULT,(void(*)(void))matQ_mult_shell); CHKERRQ(ierr);

	Mat R;
	int C_M,C_N,C_m,C_n;
	ierr = MatGetSize(proj_all_ctx.C,&C_M,&C_N); CHKERRQ(ierr);
	ierr = MatGetLocalSize(proj_all_ctx.C,&C_m,&C_n); CHKERRQ(ierr);
	if(C_n != m) SETERRQ(PETSC_COMM_SELF,1,"data inconsistency");
	if(C_N != M) SETERRQ(PETSC_COMM_SELF,1,"data inconsistency");
	ierr = MatCreateShell(PETSC_COMM_WORLD,m,C_m,M,C_M,&proj_all_ctx,&R); CHKERRQ(ierr);
	ierr = MatShellSetOperation(R,MATOP_MULT,(void(*)(void))matR_mult_shell); CHKERRQ(ierr);

	//-QTKRg
	Vec Rg;
	ierr = VecDuplicate(f,&Rg); CHKERRQ(ierr);
	ierr = MatMult(R,proj_all_ctx.g,Rg); CHKERRQ(ierr);
	PetscReal Rg_nrm2;
	ierr = VecNorm(Rg,NORM_2,&Rg_nrm2); CHKERRQ(ierr);
	PetscPrintf(PETSC_COMM_WORLD,"\t Rg_nrm2 = %6.4e\n",Rg_nrm2);
	Vec KRg;
	ierr = VecDuplicate(f,&KRg); CHKERRQ(ierr);
	ierr = MatMult(proj_all_ctx.K,Rg,KRg); CHKERRQ(ierr);
	PetscReal KRg_nrm2;
	ierr = VecNorm(KRg,NORM_2,&KRg_nrm2); CHKERRQ(ierr);
	PetscPrintf(PETSC_COMM_WORLD,"\t KRg_nrm2 = %6.4e\n",KRg_nrm2);
	ierr = VecAXPY(f,-1.,KRg); CHKERRQ(ierr);

	//+QT*(f-QTKRg)
	PetscReal f_nrm2;
	ierr = VecNorm(f,NORM_2,&f_nrm2); CHKERRQ(ierr);
	PetscPrintf(PETSC_COMM_WORLD,"\t f_nrm2 = %6.4e\n",f_nrm2);
	Vec tmp_f;
	ierr = VecDuplicate(f,&tmp_f); CHKERRQ(ierr);
	ierr = MatMult(Q,f,tmp_f); CHKERRQ(ierr);
	PetscReal QTf_nrm2;
	ierr = VecNorm(tmp_f,NORM_2,&QTf_nrm2); CHKERRQ(ierr);
	PetscPrintf(PETSC_COMM_WORLD,"\t QTf_nrm2 = %6.4e\n",QTf_nrm2);
	ierr = VecSwap(f,tmp_f); CHKERRQ(ierr);

	//CTg
	ierr = MatMultAdd(proj_all_ctx.CT,proj_all_ctx.g,f,f); CHKERRQ(ierr);
	
	ierr = VecDestroy(&KRg); CHKERRQ(ierr);
	ierr = VecDestroy(&Rg); CHKERRQ(ierr);
	ierr = VecDestroy(&tmp_f); CHKERRQ(ierr);
	ierr = MatDestroy(&Q); CHKERRQ(ierr);
	ierr = MatDestroy(&R); CHKERRQ(ierr);

    PetscFunctionReturn(0);
  }

  PetscErrorCode _postProcess_ctx_SNESSetJacobian(Mat B) {
    PetscFunctionBegin;
    ierr = MatCopy(proj_all_ctx.K,B,DIFFERENT_NONZERO_PATTERN); CHKERRQ(ierr);
    ierr = MatAXPY(B,alpha,proj_all_ctx.CTC,DIFFERENT_NONZERO_PATTERN); CHKERRQ(ierr);
    PetscFunctionReturn(0);
  }

};

struct FEMethod_DriverComplexForLazy_MaterialProjected: public  FEMethod_DriverComplexForLazy_Material,FEMethod_DriverComplexForLazy_Projected {

  FEMethod_DriverComplexForLazy_MaterialProjected(FieldInterface& _mField,matPROJ_ctx &_proj_all_ctx,BaseDirihletBC *_dirihlet_bc_method_ptr,double _lambda,double _mu,int _verbose = 0):
    FEMethod_ComplexForLazy_Data(_mField,_dirihlet_bc_method_ptr,_verbose), 
    FEMethod_DriverComplexForLazy_Material(_mField,_dirihlet_bc_method_ptr,_lambda,_mu,_verbose),
    FEMethod_DriverComplexForLazy_Projected(_mField,_proj_all_ctx,_dirihlet_bc_method_ptr,"MATERIAL_MECHANICS") {
      type_of_analysis = material_analysis;
    }

  PetscErrorCode preProcess() {
    PetscFunctionBegin;

    ierr = FEMethod_DriverComplexForLazy_Material::preProcess(); CHKERRQ(ierr);

    switch(snes_ctx) {
      case ctx_SNESSetFunction: { 
	ierr = _preProcess_ctx_SNESSetFunction(snes_x,snes_f); CHKERRQ(ierr);
      } break;
      case ctx_SNESSetJacobian: {
      } break;
      default:
	break;
    }

    switch(snes_ctx) {
      case ctx_SNESSetFunction: { 
	ierr = VecZeroEntries(snes_f); CHKERRQ(ierr);
      }
      break;
      case ctx_SNESSetJacobian:
	ierr = MatZeroEntries(proj_all_ctx.K); CHKERRQ(ierr);
	break;
      default:
	SETERRQ(PETSC_COMM_SELF,1,"not implemented");
    }

    PetscFunctionReturn(0);
  }

  PetscErrorCode operator()() {
    PetscFunctionBegin;
    ierr = FEMethod_DriverComplexForLazy_Material::operator()(&(proj_all_ctx.K)); CHKERRQ(ierr);
    PetscFunctionReturn(0);
  }

  PetscErrorCode postProcess() {
    PetscFunctionBegin;

    ierr = FEMethod_DriverComplexForLazy_Material::postProcess(); CHKERRQ(ierr);

    switch(snes_ctx) {
      case ctx_SNESSetFunction: { 
	ierr = _postProcess_ctx_SNESSetFunction(snes_f); CHKERRQ(ierr);
      }
      break;
      case ctx_SNESSetJacobian:
	ierr = MatAssemblyBegin(proj_all_ctx.K,MAT_FLUSH_ASSEMBLY); CHKERRQ(ierr);
	ierr = MatAssemblyEnd(proj_all_ctx.K,MAT_FLUSH_ASSEMBLY); CHKERRQ(ierr);
	ierr = dirihlet_bc_method_ptr->SetDirihletBC_to_MatrixDiagonal(this,proj_all_ctx.K); CHKERRQ(ierr);
	ierr = MatAssemblyBegin(proj_all_ctx.K,MAT_FINAL_ASSEMBLY); CHKERRQ(ierr);
	ierr = MatAssemblyEnd(proj_all_ctx.K,MAT_FINAL_ASSEMBLY); CHKERRQ(ierr);
	ierr = _postProcess_ctx_SNESSetJacobian(*snes_B); CHKERRQ(ierr);
	break;
      default:
	SETERRQ(PETSC_COMM_SELF,1,"not implemented");
    }

    PetscFunctionReturn(0);
  }

};

struct FEMethod_DriverComplexForLazy_MeshSmoothingProjected: public FEMethod_DriverComplexForLazy_MeshSmoothing,FEMethod_DriverComplexForLazy_Projected {

  FEMethod_DriverComplexForLazy_MeshSmoothingProjected(FieldInterface& _mField,matPROJ_ctx &_proj_all_ctx,BaseDirihletBC *_dirihlet_bc_method_ptr,int _verbose = 0):
    FEMethod_ComplexForLazy_Data(_mField,_dirihlet_bc_method_ptr,_verbose), 
    FEMethod_DriverComplexForLazy_MeshSmoothing(_mField,_dirihlet_bc_method_ptr,_verbose),
    FEMethod_DriverComplexForLazy_Projected(_mField,_proj_all_ctx,_dirihlet_bc_method_ptr,"MESH_SMOOTHING") {}

  PetscErrorCode preProcess() {
    PetscFunctionBegin;

    ierr = FEMethod_DriverComplexForLazy_MeshSmoothing::preProcess(); CHKERRQ(ierr);

    switch(snes_ctx) {
      case ctx_SNESSetFunction: { 
	ierr = _preProcess_ctx_SNESSetFunction(snes_x,snes_f); CHKERRQ(ierr);
      } break;
      case ctx_SNESSetJacobian: {
      } break;
      default:
	break;
    }

    switch(snes_ctx) {
      case ctx_SNESSetFunction: { 
	ierr = VecZeroEntries(snes_f); CHKERRQ(ierr);
      }
      break;
      case ctx_SNESSetJacobian:
	ierr = MatZeroEntries(proj_all_ctx.K); CHKERRQ(ierr);
	break;
      default:
	SETERRQ(PETSC_COMM_SELF,1,"not implemented");
    }

    PetscFunctionReturn(0);
  }

  PetscErrorCode operator()() {
    PetscFunctionBegin;
    ierr = FEMethod_DriverComplexForLazy_MeshSmoothing::operator()(&(proj_all_ctx.K)); CHKERRQ(ierr);
    PetscFunctionReturn(0);
  }

  PetscErrorCode postProcess() {
    PetscFunctionBegin;

    ierr = FEMethod_DriverComplexForLazy_Material::postProcess(); CHKERRQ(ierr);

    switch(snes_ctx) {
      case ctx_SNESSetFunction: { 
	ierr = _postProcess_ctx_SNESSetFunction(snes_f); CHKERRQ(ierr);
      }
      break;
      case ctx_SNESSetJacobian:
	ierr = MatAssemblyBegin(proj_all_ctx.K,MAT_FLUSH_ASSEMBLY); CHKERRQ(ierr);
	ierr = MatAssemblyEnd(proj_all_ctx.K,MAT_FLUSH_ASSEMBLY); CHKERRQ(ierr);
	ierr = dirihlet_bc_method_ptr->SetDirihletBC_to_MatrixDiagonal(this,proj_all_ctx.K); CHKERRQ(ierr);
	ierr = MatAssemblyBegin(proj_all_ctx.K,MAT_FINAL_ASSEMBLY); CHKERRQ(ierr);
	ierr = MatAssemblyEnd(proj_all_ctx.K,MAT_FINAL_ASSEMBLY); CHKERRQ(ierr);
	ierr = _postProcess_ctx_SNESSetJacobian(*snes_B); CHKERRQ(ierr);
	break;
      default:
	SETERRQ(PETSC_COMM_SELF,1,"not implemented");
    }

    PetscFunctionReturn(0);
  }

};

struct FEMethod_DriverComplexForLazy_CoupledSpatial: public FEMethod_DriverComplexForLazy_Spatial {

  matPROJ_ctx& proj_all_ctx;
  FEMethod_DriverComplexForLazy_CoupledSpatial(FieldInterface& _mField,matPROJ_ctx &_proj_all_ctx,BaseDirihletBC *_dirihlet_bc_method_ptr,double _lambda,double _mu,int _verbose = 0):
    FEMethod_ComplexForLazy_Data(_mField,_dirihlet_bc_method_ptr,_verbose), 
    FEMethod_DriverComplexForLazy_Spatial(_mField,_dirihlet_bc_method_ptr,_lambda,_mu,_verbose),
    proj_all_ctx(_proj_all_ctx) {}


  PetscErrorCode preProcess() {
    PetscFunctionBegin;

    //PetscSynchronizedPrintf(PETSC_COMM_WORLD,"Start Assembly\n");
    ierr = PetscGetTime(&v1); CHKERRQ(ierr);
    ierr = PetscGetCPUTime(&t1); CHKERRQ(ierr);
    switch(snes_ctx) {
      case ctx_SNESNone:
      case ctx_SNESSetFunction: { 
	ierr = VecAssemblyBegin(snes_f); CHKERRQ(ierr);
	ierr = VecAssemblyEnd(snes_f); CHKERRQ(ierr);
      }
      break;
      case ctx_SNESSetJacobian: {
	ierr = MatAssemblyBegin(proj_all_ctx.K,MAT_FLUSH_ASSEMBLY); CHKERRQ(ierr);
	ierr = MatAssemblyEnd(proj_all_ctx.K,MAT_FLUSH_ASSEMBLY); CHKERRQ(ierr);
      }
      break;
      default:
	SETERRQ(PETSC_COMM_SELF,1,"not implemented");
    }

    PetscFunctionReturn(0);
  }

  virtual PetscErrorCode CalculateSpatialCoupledTangent(Mat B) {
    PetscFunctionBegin;

    switch(snes_ctx) {
      case ctx_SNESNone:
      case ctx_SNESSetJacobian:
	ierr = MatSetValues(B,
	  RowGlobMaterial[i_nodes].size(),&*(RowGlobMaterial[i_nodes].begin()),
	  ColGlobSpatial[i_nodes].size(),&*(ColGlobSpatial[i_nodes].begin()),
	  &*(KHh.data().begin()),ADD_VALUES); CHKERRQ(ierr);
	for(int ee = 0;ee<6;ee++) {
	  ierr = MatSetValues(B,
	    RowGlobMaterial[i_nodes].size(),&*(RowGlobMaterial[i_nodes].begin()),
	    ColGlobSpatial[1+ee].size(),&*(ColGlobSpatial[1+ee].begin()),
	    &*(KHedge_data[ee].data().begin()),ADD_VALUES); CHKERRQ(ierr);
	}
	for(int ff = 0;ff<4;ff++) {
	  ierr = MatSetValues(B,
	    RowGlobMaterial[i_nodes].size(),&*(RowGlobMaterial[i_nodes].begin()),
	    ColGlobSpatial[1+6+ff].size(),&*(ColGlobSpatial[1+6+ff].begin()),
	    &*(KHface_data[ff].data().begin()),ADD_VALUES); CHKERRQ(ierr);
	}
	ierr = MatSetValues(B,
	  RowGlobMaterial[i_nodes].size(),&*(RowGlobMaterial[i_nodes].begin()),
	  ColGlobSpatial[i_volume].size(),&*(ColGlobSpatial[i_volume].begin()),
	  &*(KHvolume.data().begin()),ADD_VALUES); CHKERRQ(ierr);
	break;
      default:
	SETERRQ(PETSC_COMM_SELF,1,"not implemented");
    }

    PetscFunctionReturn(0);

  }

  PetscErrorCode operator()() {
    PetscFunctionBegin;

    ierr = OpComplexForLazyStart(); CHKERRQ(ierr);
    ierr = GetIndicesSpatial(); CHKERRQ(ierr);

    ierr = dirihlet_bc_method_ptr->SetDirihletBC_to_ElementIndicies(this,RowGlobSpatial,ColGlobSpatial,DirihletBC); CHKERRQ(ierr);

    switch(snes_ctx) {
      case ctx_SNESSetFunction: { 
	ierr = CalculateSpatialFint(snes_f); CHKERRQ(ierr);
      }
      break;
      case ctx_SNESSetJacobian:
	ierr = CalculateSpatialTangent(proj_all_ctx.K); CHKERRQ(ierr);
	ierr = GetIndicesRow(RowGlobMaterial,material_field_name); CHKERRQ(ierr);
	ierr = dirihlet_bc_method_ptr->SetDirihletBC_to_ElementIndiciesRow(this,RowGlobMaterial,DirihletBC); CHKERRQ(ierr);
	ierr = CalculateSpatialCoupledTangent(proj_all_ctx.K); CHKERRQ(ierr);
	break;
      default:
	SETERRQ(PETSC_COMM_SELF,1,"not implemented");
    }

    for(_IT_CUBITMESHSETS_BY_BCDATA_TYPE_FOR_LOOP_(mField,SideSet|PressureSet,it)) {

      pressure_cubit_bc_data mydata;
      ierr = it->get_cubit_bc_data_structure(mydata); CHKERRQ(ierr);

      double t_val_ = *(this->t_val)*mydata.data.value1;
      double t[] = { 0,0,-t_val_, 0,0,-t_val_, 0,0,-t_val_ };

      Range NeumannSideSet;
      ierr = it->get_Cubit_msId_entities_by_dimension(mField.get_moab(),2,NeumannSideSet,true); CHKERRQ(ierr);

      switch(snes_ctx) {
	case ctx_SNESSetFunction: { 
	  ierr = CaluclateSpatialFext(snes_f,t,NeumannSideSet); CHKERRQ(ierr);
	}
	break;
	case ctx_SNESSetJacobian:
	  ierr = CalculateSpatialTangentExt(proj_all_ctx.K,t,NeumannSideSet); CHKERRQ(ierr);
	  break;
	default:
	  SETERRQ(PETSC_COMM_SELF,1,"not implemented");
      }

    }

    PetscFunctionReturn(0);
  }

  PetscErrorCode postProcess() {
    PetscFunctionBegin;

    switch(snes_ctx) {
      case ctx_SNESSetFunction: { 
	ierr = VecAssemblyBegin(snes_f); CHKERRQ(ierr);
	ierr = VecAssemblyEnd(snes_f); CHKERRQ(ierr);
	ierr = dirihlet_bc_method_ptr->SetDirihletBC_to_RHS(this,snes_f); CHKERRQ(ierr);
      }
      break;
      case ctx_SNESSetJacobian: {
	ierr = MatAssemblyBegin(proj_all_ctx.K,MAT_FLUSH_ASSEMBLY); CHKERRQ(ierr);
	ierr = MatAssemblyEnd(proj_all_ctx.K,MAT_FLUSH_ASSEMBLY); CHKERRQ(ierr);
	ierr = dirihlet_bc_method_ptr->SetDirihletBC_to_MatrixDiagonal(this,proj_all_ctx.K); CHKERRQ(ierr);
      }
      break;
      default:
	SETERRQ(PETSC_COMM_SELF,1,"not implemented");
    }

    PetscFunctionReturn(0);
  }


};

struct FEMethod_DriverComplexForLazy_CoupledMaterial: public FEMethod_DriverComplexForLazy_Material {

  matPROJ_ctx& proj_all_ctx;
  FEMethod_DriverComplexForLazy_CoupledMaterial(FieldInterface& _mField,matPROJ_ctx &_proj_all_ctx,BaseDirihletBC *_dirihlet_bc_method_ptr,double _lambda,double _mu,int _verbose = 0):
    FEMethod_ComplexForLazy_Data(_mField,_dirihlet_bc_method_ptr,_verbose), 
    FEMethod_DriverComplexForLazy_Material(_mField,_dirihlet_bc_method_ptr,_lambda,_mu,_verbose),
    proj_all_ctx(_proj_all_ctx) {}

  PetscErrorCode preProcess() {
    PetscFunctionBegin;

    //PetscSynchronizedPrintf(PETSC_COMM_WORLD,"Start Assembly\n");
    ierr = PetscGetTime(&v1); CHKERRQ(ierr);
    ierr = PetscGetCPUTime(&t1); CHKERRQ(ierr);
    switch(snes_ctx) {
      case ctx_SNESNone:
      case ctx_SNESSetFunction: { 
	ierr = VecAssemblyBegin(snes_f); CHKERRQ(ierr);
	ierr = VecAssemblyEnd(snes_f); CHKERRQ(ierr);
      }
      break;
      case ctx_SNESSetJacobian: {
	ierr = MatAssemblyBegin(proj_all_ctx.K,MAT_FLUSH_ASSEMBLY); CHKERRQ(ierr);
	ierr = MatAssemblyEnd(proj_all_ctx.K,MAT_FLUSH_ASSEMBLY); CHKERRQ(ierr);
      }
      break;
      default:
	SETERRQ(PETSC_COMM_SELF,1,"not implemented");
    }

    PetscFunctionReturn(0);
  }

  virtual PetscErrorCode CalculateMaterialCoupledTangent(Mat B) {
    PetscFunctionBegin;

    switch(snes_ctx) {
      case ctx_SNESNone:
      case ctx_SNESSetJacobian:
	ierr = MatSetValues(B,
	  RowGlobSpatial[i_nodes].size(),&*(RowGlobSpatial[i_nodes].begin()),
	  ColGlobMaterial[i_nodes].size(),&*(ColGlobMaterial[i_nodes].begin()),
	  &*(KhH.data().begin()),ADD_VALUES); CHKERRQ(ierr);
	for(int ee = 0;ee<6;ee++) {
	  ierr = MatSetValues(B,
	    RowGlobSpatial[i_nodes].size(),&*(RowGlobSpatial[i_nodes].begin()),
	    ColGlobMaterial[1+ee].size(),&*(ColGlobMaterial[1+ee].begin()),
	    &*(KedgeH_data[ee].data().begin()),ADD_VALUES); CHKERRQ(ierr);
	}
	for(int ff = 0;ff<4;ff++) {
	  ierr = MatSetValues(B,
	    RowGlobSpatial[i_nodes].size(),&*(RowGlobSpatial[i_nodes].begin()),
	    ColGlobMaterial[1+6+ff].size(),&*(ColGlobMaterial[1+6+ff].begin()),
	    &*(KfaceH_data[ff].data().begin()),ADD_VALUES); CHKERRQ(ierr);
	}
	ierr = MatSetValues(B,
	  RowGlobSpatial[i_nodes].size(),&*(RowGlobSpatial[i_nodes].begin()),
	  ColGlobMaterial[i_volume].size(),&*(ColGlobMaterial[i_volume].begin()),
	  &*(KvolumeH.data().begin()),ADD_VALUES); CHKERRQ(ierr);
	break;
      default:
	SETERRQ(PETSC_COMM_SELF,1,"not implemented");
    }

    PetscFunctionReturn(0);

  }

  PetscErrorCode operator()() {
    PetscFunctionBegin;
    ierr = FEMethod_DriverComplexForLazy_Material::operator()(&(proj_all_ctx.K)); CHKERRQ(ierr);

    switch(snes_ctx) {
      case ctx_SNESSetFunction:  
      break;
      case ctx_SNESSetJacobian:
	ierr = GetIndicesRow(RowGlobSpatial,spatial_field_name); CHKERRQ(ierr);
	ierr = dirihlet_bc_method_ptr->SetDirihletBC_to_ElementIndiciesRow(this,RowGlobSpatial,DirihletBC); CHKERRQ(ierr);
	ierr = CalculateMaterialCoupledTangent(proj_all_ctx.K); CHKERRQ(ierr);
	break;
      default:
	SETERRQ(PETSC_COMM_SELF,1,"not implemented");
    }

    PetscFunctionReturn(0);
  }

  PetscErrorCode postProcess() {
    PetscFunctionBegin;

    switch(snes_ctx) {
      case ctx_SNESSetFunction: { 
	ierr = VecAssemblyBegin(snes_f); CHKERRQ(ierr);
	ierr = VecAssemblyEnd(snes_f); CHKERRQ(ierr);
	ierr = dirihlet_bc_method_ptr->SetDirihletBC_to_RHS(this,snes_f); CHKERRQ(ierr);
      }
      break;
      case ctx_SNESSetJacobian: {
	ierr = MatAssemblyBegin(proj_all_ctx.K,MAT_FLUSH_ASSEMBLY); CHKERRQ(ierr);
	ierr = MatAssemblyEnd(proj_all_ctx.K,MAT_FLUSH_ASSEMBLY); CHKERRQ(ierr);
	ierr = dirihlet_bc_method_ptr->SetDirihletBC_to_MatrixDiagonal(this,proj_all_ctx.K); CHKERRQ(ierr);
      }
      break;
      default:
	SETERRQ(PETSC_COMM_SELF,1,"not implemented");
    }

    PetscFunctionReturn(0);
  }

};

struct FEMethod_DriverComplexForLazy_CoupledMeshSmoother: public FEMethod_DriverComplexForLazy_MeshSmoothing {

  matPROJ_ctx& proj_all_ctx;
  FEMethod_DriverComplexForLazy_CoupledMeshSmoother(FieldInterface& _mField,matPROJ_ctx &_proj_all_ctx,BaseDirihletBC *_dirihlet_bc_method_ptr,int _verbose = 0):
    FEMethod_ComplexForLazy_Data(_mField,_dirihlet_bc_method_ptr,_verbose), 
    FEMethod_DriverComplexForLazy_MeshSmoothing(_mField,_dirihlet_bc_method_ptr,_verbose),
    proj_all_ctx(_proj_all_ctx) {}

  PetscErrorCode preProcess() {
    PetscFunctionBegin;

    //PetscSynchronizedPrintf(PETSC_COMM_WORLD,"Start Assembly\n");
    ierr = PetscGetTime(&v1); CHKERRQ(ierr);
    ierr = PetscGetCPUTime(&t1); CHKERRQ(ierr);
    switch(snes_ctx) {
      case ctx_SNESNone:
      case ctx_SNESSetFunction: { 
	ierr = VecAssemblyBegin(snes_f); CHKERRQ(ierr);
	ierr = VecAssemblyEnd(snes_f); CHKERRQ(ierr);
      }
      break;
      case ctx_SNESSetJacobian: {
	ierr = MatAssemblyBegin(proj_all_ctx.K,MAT_FLUSH_ASSEMBLY); CHKERRQ(ierr);
	ierr = MatAssemblyEnd(proj_all_ctx.K,MAT_FLUSH_ASSEMBLY); CHKERRQ(ierr);
      }
      break;
      default:
	SETERRQ(PETSC_COMM_SELF,1,"not implemented");
    }

    PetscFunctionReturn(0);
  }

  PetscErrorCode operator()() {
    PetscFunctionBegin;
    ierr = FEMethod_DriverComplexForLazy_MeshSmoothing::operator()(&(proj_all_ctx.K)); CHKERRQ(ierr);
    PetscFunctionReturn(0);
  }

  PetscErrorCode postProcess() {
    PetscFunctionBegin;

    switch(snes_ctx) {
      case ctx_SNESSetFunction: { 
	ierr = VecAssemblyBegin(snes_f); CHKERRQ(ierr);
	ierr = VecAssemblyEnd(snes_f); CHKERRQ(ierr);
	ierr = dirihlet_bc_method_ptr->SetDirihletBC_to_RHS(this,snes_f); CHKERRQ(ierr);
      }
      break;
      case ctx_SNESSetJacobian: {
	ierr = MatAssemblyBegin(proj_all_ctx.K,MAT_FLUSH_ASSEMBLY); CHKERRQ(ierr);
	ierr = MatAssemblyEnd(proj_all_ctx.K,MAT_FLUSH_ASSEMBLY); CHKERRQ(ierr);
	ierr = dirihlet_bc_method_ptr->SetDirihletBC_to_MatrixDiagonal(this,proj_all_ctx.K); CHKERRQ(ierr);
      }
      break;
      default:
	SETERRQ(PETSC_COMM_SELF,1,"not implemented");
    }

    PetscFunctionReturn(0);
  }

};

struct FEMethod_DriverComplexForLazy_CoupledProjected: public  FEMethod_DriverComplexForLazy_Projected {

  FEMethod_DriverComplexForLazy_CoupledProjected(FieldInterface& _mField,matPROJ_ctx &_proj_all_ctx,BaseDirihletBC *_dirihlet_bc_method_ptr,string _problem_name,int _verbose = 0):
    FEMethod_ComplexForLazy_Data(_mField,_dirihlet_bc_method_ptr,_verbose), 
    FEMethod_DriverComplexForLazy_Projected(_mField,_proj_all_ctx,_dirihlet_bc_method_ptr,_problem_name) {} 

  PetscErrorCode preProcess() {
    PetscFunctionBegin;

    switch(snes_ctx) {
      case ctx_SNESSetFunction: { 
	ierr = _preProcess_ctx_SNESSetFunction(snes_x,snes_f); CHKERRQ(ierr);
      } break;
      case ctx_SNESSetJacobian: {
      } break;
      default:
	break;
    }

    switch(snes_ctx) {
      case ctx_SNESSetFunction: { 
	ierr = VecZeroEntries(snes_f); CHKERRQ(ierr);
      }
      break;
      case ctx_SNESSetJacobian:
	ierr = MatZeroEntries(proj_all_ctx.K); CHKERRQ(ierr);
	break;
      default:
	SETERRQ(PETSC_COMM_SELF,1,"not implemented");
    }


    PetscFunctionReturn(0);
  }

  PetscErrorCode postProcess() {
    PetscFunctionBegin;

    switch(snes_ctx) {
      case ctx_SNESSetFunction: { 
	ierr = VecAssemblyBegin(snes_f); CHKERRQ(ierr);
	ierr = VecAssemblyEnd(snes_f); CHKERRQ(ierr);
	ierr = _postProcess_ctx_SNESSetFunction(snes_f); CHKERRQ(ierr);
	ierr = VecAssemblyBegin(snes_f); CHKERRQ(ierr);
	ierr = VecAssemblyEnd(snes_f); CHKERRQ(ierr);
      }
      break;
      case ctx_SNESSetJacobian:
	ierr = MatAssemblyBegin(proj_all_ctx.K,MAT_FINAL_ASSEMBLY); CHKERRQ(ierr);
	ierr = MatAssemblyEnd(proj_all_ctx.K,MAT_FINAL_ASSEMBLY); CHKERRQ(ierr);
	ierr = _postProcess_ctx_SNESSetJacobian(*snes_B); CHKERRQ(ierr);
	break;
      default:
	SETERRQ(PETSC_COMM_SELF,1,"not implemented");
    }

    PetscFunctionReturn(0);
  }

};

}

#endif

