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

#include "SnesCtx.hpp"
#include "FEMethod_ComplexForLazy.hpp"
#include "ArcLengthTools.hpp"

#include "petscShellMATs_ConstrainsByMarkAinsworth.hpp"
#include "FEMethod_SurfaceConstrains.hpp"

extern "C" {
#include "gm_rule.h"
}

using namespace boost::numeric;

namespace MoFEM {

struct NonLinearSpatialElasticFEMthod: public FEMethod_ComplexForLazy {

  ArcLengthCtx* arcPtr;
  bool isCoupledProblem;

  NonLinearSpatialElasticFEMthod(FieldInterface& _mField,double _lambda,double _mu,int _verbose = 0): 
    FEMethod_ComplexForLazy_Data(_mField,_verbose),
    FEMethod_ComplexForLazy(_mField,FEMethod_ComplexForLazy::spatail_analysis,_lambda,_mu,0,_verbose),arcPtr(NULL),
    isCoupledProblem(false) {}

  double *thermalLoadFactor;
  //set load factor
  PetscErrorCode set_thermal_load_factor(double t_thermal_load_factor_val_) {
      PetscFunctionBegin;
      *thermalLoadFactor = t_thermal_load_factor_val_;
      PetscFunctionReturn(0);
  }
  Tag thThermalLoadFactor;

  NonLinearSpatialElasticFEMthod(FieldInterface& _mField,double _lambda,double _mu,ArcLengthCtx* arc_ptr,int _verbose = 0): 
    FEMethod_ComplexForLazy_Data(_mField,_verbose),
    FEMethod_ComplexForLazy(_mField,FEMethod_ComplexForLazy::spatail_analysis,_lambda,_mu,0,_verbose),arcPtr(arc_ptr),
    isCoupledProblem(false) {

    double def_t_val = 0;
    const EntityHandle root_meshset = mField.get_moab().get_root_set();
    rval = mField.get_moab().tag_get_handle("_ThermalExpansionFactor_t_alpha_val",1,MB_TYPE_DOUBLE,thThermalLoadFactor,MB_TAG_CREAT|MB_TAG_EXCL|MB_TAG_MESH,&def_t_val); 
    if(rval == MB_ALREADY_ALLOCATED) {
      rval = mField.get_moab().tag_get_by_ptr(thThermalLoadFactor,&root_meshset,1,(const void**)&thermalLoadFactor); CHKERR_THROW(rval);
    } else {
      CHKERR_THROW(rval);
      rval = mField.get_moab().tag_set_data(thThermalLoadFactor,&root_meshset,1,&def_t_val); CHKERR_THROW(rval);
      rval = mField.get_moab().tag_get_by_ptr(thThermalLoadFactor,&root_meshset,1,(const void**)&thermalLoadFactor); CHKERR_THROW(rval);
    }

  }

  virtual PetscErrorCode AssembleSpatialFint(Vec f) {
    PetscFunctionBegin;

    switch(snes_ctx) {
      case CTX_SNESSETFUNCTION: { 
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
      break;
    }

    if(arcPtr != NULL) {

      switch(snes_ctx) {
        case CTX_SNESNONE:
        case CTX_SNESSETFUNCTION: { 
  	  analysis _type_of_analysis = type_of_analysis;
  	  type_of_analysis = scaled_themp_direvative_spatial;
  	  ierr = GetFint(); CHKERRQ(ierr);
  	  iFint_h /= -eps;
  	  ierr = VecSetValues(arcPtr->F_lambda,RowGlobSpatial[i_nodes].size(),&(RowGlobSpatial[i_nodes][0]),&(iFint_h.data()[0]),ADD_VALUES); CHKERRQ(ierr);
  	  for(int ee = 0;ee<6;ee++) {
  	    if(RowGlobSpatial[1+ee].size()>0) {
  	      iFint_h_edge_data[ee] /= -eps;
  	      ierr = VecSetValues(arcPtr->F_lambda,RowGlobSpatial[1+ee].size(),&(RowGlobSpatial[1+ee][0]),&(iFint_h_edge_data[ee].data()[0]),ADD_VALUES); CHKERRQ(ierr);
  	    }
  	  }
  	  for(int ff = 0;ff<4;ff++) {
  	    if(RowGlobSpatial[1+6+ff].size()>0) {
  	      iFint_h_face_data[ff] /= -eps;
  	      ierr = VecSetValues(arcPtr->F_lambda,RowGlobSpatial[1+6+ff].size(),&(RowGlobSpatial[1+6+ff][0]),&(iFint_h_face_data[ff].data()[0]),ADD_VALUES); CHKERRQ(ierr);
  	    }
  	  }
  	  if(RowGlobSpatial[i_volume].size()>0) {
  	    iFint_h_volume /= -eps;
  	    ierr = VecSetValues(arcPtr->F_lambda,RowGlobSpatial[i_volume].size(),&(RowGlobSpatial[i_volume][0]),&(iFint_h_volume.data()[0]),ADD_VALUES); CHKERRQ(ierr);
  	  }
  	  type_of_analysis = _type_of_analysis;
  	}
	break;
        default:
        break;
      }

    }

    PetscFunctionReturn(0);
  }

  virtual PetscErrorCode AssembleSpatialTangent(Mat B) {
    PetscFunctionBegin;

    switch(snes_ctx) {
      case CTX_SNESNONE:
      case CTX_SNESSETFUNCTION:
      case CTX_SNESSETJACOBIAN:
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

  virtual PetscErrorCode AssembleSpatialCoupledTangent(Mat B) {
    PetscFunctionBegin;
    switch(snes_ctx) {
      case CTX_SNESSETJACOBIAN:
	if(KHh.size1()!=RowGlobMaterial[0].size()) {
	  SETERRQ(PETSC_COMM_SELF,1,"KHh.size()!=RowGlobMaterial[0].size()");
	}
        ierr = MatSetValues(B,
	  RowGlobMaterial[0].size(),&*(RowGlobMaterial[0].begin()),
	  ColGlobSpatial[i_nodes].size(),&*(ColGlobSpatial[i_nodes].begin()),
	  &*(KHh.data().begin()),ADD_VALUES); CHKERRQ(ierr);
        for(int ee = 0;ee<6;ee++) {
	  ierr = MatSetValues(B,
	    RowGlobMaterial[0].size(),&*(RowGlobMaterial[0].begin()),
	    ColGlobSpatial[1+ee].size(),&*(ColGlobSpatial[1+ee].begin()),
	    &*(KHedge_data[ee].data().begin()),ADD_VALUES); CHKERRQ(ierr);
        }
        for(int ff = 0;ff<4;ff++) {
	  ierr = MatSetValues(B,
	    RowGlobMaterial[0].size(),&*(RowGlobMaterial[0].begin()),
	    ColGlobSpatial[1+6+ff].size(),&*(ColGlobSpatial[1+6+ff].begin()),
	    &*(KHface_data[ff].data().begin()),ADD_VALUES); CHKERRQ(ierr);
        }
        ierr = MatSetValues(B,
	  RowGlobMaterial[0].size(),&*(RowGlobMaterial[0].begin()),
	  ColGlobSpatial[i_volume].size(),&*(ColGlobSpatial[i_volume].begin()),
	  &*(KHvolume.data().begin()),ADD_VALUES); CHKERRQ(ierr);
        break;
      default:
        SETERRQ(PETSC_COMM_SELF,1,"not implemented");
    }
    PetscFunctionReturn(0);
  }

  PetscErrorCode preProcess() {
    PetscFunctionBegin;
    if(arcPtr!=NULL) {
      ierr = set_thermal_load_factor(arcPtr->get_FieldData()); CHKERRQ(ierr);
      thermal_load_factor = *thermalLoadFactor;
    }
    PetscFunctionReturn(0);
  }

  ublas::matrix<double> gaussPts;
  PetscErrorCode operator()() {
    PetscFunctionBegin;
    //gauss points
    int order = 1;
    for(_IT_GET_FEROW_BY_TYPE_DOFS_FOR_LOOP_(this,spatial_field_name,MBEDGE,dit)) {
      order = max(order,dit->get_max_order());
    }
    int order_themp = 0;
    for(_IT_GET_FEDATA_BY_TYPE_DOFS_FOR_LOOP_(this,termal_field_name,MBVERTEX,dit)) {
      order_themp = 1;
      break;
    }
    for(_IT_GET_FEDATA_BY_TYPE_DOFS_FOR_LOOP_(this,termal_field_name,MBEDGE,dit)) {
      order_themp = max(order_themp,dit->get_max_order());
    }
    int rule = order-1 + order_themp > 0 ? order-1 + order_themp : 0; //FIXME rule with temperature
    int nb_gauss_pts = gm_rule_size(rule,3);
    gaussPts.resize(4,nb_gauss_pts);
    ierr = Grundmann_Moeller_integration_points_3D_TET(
      rule,&gaussPts(0,0),&gaussPts(1,0),&gaussPts(2,0),&gaussPts(3,0)); CHKERRQ(ierr);
    g_NTET.resize(4*nb_gauss_pts);
    ierr = ShapeMBTET(&g_NTET[0],&gaussPts(0,0),&gaussPts(1,0),&gaussPts(2,0),nb_gauss_pts); CHKERRQ(ierr);
    g_TET_W = &(gaussPts(3,0));
    ierr = OpComplexForLazyStart(); CHKERRQ(ierr);
    ierr = GetIndicesSpatial(); CHKERRQ(ierr);
    switch(snes_ctx) {
      case CTX_SNESNONE:
      case CTX_SNESSETFUNCTION: { 
        ierr = GetFint(); CHKERRQ(ierr);
	ierr = AssembleSpatialFint(snes_f); CHKERRQ(ierr);
      }
      break;
      case CTX_SNESSETJACOBIAN:
	ierr = GetTangent(); CHKERRQ(ierr);
	ierr = AssembleSpatialTangent(*snes_B); CHKERRQ(ierr);
	if(isCoupledProblem) {
	  ierr = GetIndicesRow(RowGlobMaterial,material_field_name); CHKERRQ(ierr);
	  ierr = AssembleSpatialCoupledTangent(*snes_B); CHKERRQ(ierr);
	}
	break;
      default:
	SETERRQ(PETSC_COMM_SELF,1,"not implemented");
    }
    PetscFunctionReturn(0);
  }

  PetscErrorCode postProcess() {
    PetscFunctionBegin;
    switch(snes_ctx) {
      case CTX_SNESSETFUNCTION: { 
	ierr = VecAssemblyBegin(snes_f); CHKERRQ(ierr);
	ierr = VecAssemblyEnd(snes_f); CHKERRQ(ierr);
      }
      case CTX_SNESNONE: {
	if(arcPtr!=NULL) {
	  ierr = VecAssemblyBegin(arcPtr->F_lambda); CHKERRQ(ierr);
	  ierr = VecAssemblyEnd(arcPtr->F_lambda); CHKERRQ(ierr);
	}
      }
      break;
      case CTX_SNESSETJACOBIAN: {
	ierr = MatAssemblyBegin(*snes_B,MAT_FLUSH_ASSEMBLY); CHKERRQ(ierr);
	ierr = MatAssemblyEnd(*snes_B,MAT_FLUSH_ASSEMBLY); CHKERRQ(ierr);
      }
      break;
      default:
	SETERRQ(PETSC_COMM_SELF,1,"not implemented");
    }
    PetscFunctionReturn(0);
  }

};

struct EshelbyFEMethod: public NonLinearSpatialElasticFEMthod {

  EshelbyFEMethod(FieldInterface& _mField,double _lambda,double _mu,int _verbose = 0):
    FEMethod_ComplexForLazy_Data(_mField,_verbose), 
    NonLinearSpatialElasticFEMthod(_mField,_lambda,_mu,_verbose) {
    type_of_analysis = material_analysis;
  }

  EshelbyFEMethod(FieldInterface& _mField,double _lambda,double _mu,ArcLengthCtx* _arc_ptr,int _verbose = 0):
    FEMethod_ComplexForLazy_Data(_mField,_verbose), NonLinearSpatialElasticFEMthod(_mField,_lambda,_mu,_arc_ptr,_verbose) {
    type_of_analysis = material_analysis;
  } 

  virtual PetscErrorCode AssembleMaterialTangent(Mat B) {
    PetscFunctionBegin;
    switch(snes_ctx) {
      case CTX_SNESNONE:
      case CTX_SNESSETFUNCTION:
      case CTX_SNESSETJACOBIAN:
	ierr = MatSetValues(B,
	  RowGlobMaterial[0].size(),&*(RowGlobMaterial[0].begin()),
	  ColGlobMaterial[0].size(),&*(ColGlobMaterial[0].begin()),
	  &*(KHH.data().begin()),ADD_VALUES); CHKERRQ(ierr);
	break;
      default:
	SETERRQ(PETSC_COMM_SELF,1,"not implemented");
    }
    PetscFunctionReturn(0);
  }

  virtual PetscErrorCode AssembleMaterialCoupledTangent(Mat B) {
    PetscFunctionBegin;
    switch(snes_ctx) {
      case CTX_SNESSETJACOBIAN:
	if(RowGlobSpatial.empty()) {
	  SETERRQ(PETSC_COMM_SELF,1,"RowGlobSpatial.empty()");
	}
	if(KhH.size1() != RowGlobSpatial[i_nodes].size()) {
	  SETERRQ(PETSC_COMM_SELF,1,"KhH.size1() != RowGlobSpatial[i_nodes].size()");
	}
	ierr = MatSetValues(B,
	  RowGlobSpatial[i_nodes].size(),&*(RowGlobSpatial[i_nodes].begin()),
	  ColGlobMaterial[i_nodes].size(),&*(ColGlobMaterial[i_nodes].begin()),
	  &*(KhH.data().begin()),ADD_VALUES); CHKERRQ(ierr);
	for(int ee = 0;ee<6;ee++) {
	  ierr = MatSetValues(B,
	    RowGlobSpatial[1+ee].size(),&*(RowGlobSpatial[1+ee].begin()),
	    ColGlobMaterial[i_nodes].size(),&*(ColGlobMaterial[i_nodes].begin()),
	    &*(KedgeH_data[ee].data().begin()),ADD_VALUES); CHKERRQ(ierr);
	}
	for(int ff = 0;ff<4;ff++) {
	  ierr = MatSetValues(B,
	    RowGlobSpatial[1+6+ff].size(),&*(RowGlobSpatial[1+6+ff].begin()),
	    ColGlobMaterial[i_nodes].size(),&*(ColGlobMaterial[i_nodes].begin()),
	    &*(KfaceH_data[ff].data().begin()),ADD_VALUES); CHKERRQ(ierr);
	}
	ierr = MatSetValues(B,
	  RowGlobSpatial[i_volume].size(),&*(RowGlobSpatial[i_volume].begin()),
	  ColGlobMaterial[i_nodes].size(),&*(ColGlobMaterial[i_nodes].begin()),
	  &*(KvolumeH.data().begin()),ADD_VALUES); CHKERRQ(ierr);
	break;
      default:
	SETERRQ(PETSC_COMM_SELF,1,"not implemented");
    }
    PetscFunctionReturn(0);
  }

  virtual PetscErrorCode AssembleMaterialFint(Vec f) {
    PetscFunctionBegin;
    switch(snes_ctx) {
      case CTX_SNESNONE:
      case CTX_SNESSETFUNCTION: { 
	ierr = VecSetOption(f,VEC_IGNORE_NEGATIVE_INDICES,PETSC_TRUE);  CHKERRQ(ierr);
	ierr = VecSetValues(f,RowGlobMaterial[0].size(),&(RowGlobMaterial[0][0]),&(Fint_H.data()[0]),ADD_VALUES); CHKERRQ(ierr);
      }
      break;
      default:
	SETERRQ(PETSC_COMM_SELF,1,"not implemented");
    }
    PetscFunctionReturn(0);
  }

  PetscErrorCode operator()() {
    PetscFunctionBegin;
    //gauss points
    int order = 1;
    for(_IT_GET_FEROW_BY_TYPE_DOFS_FOR_LOOP_(this,spatial_field_name,MBEDGE,dit)) {
      order = max(order,dit->get_max_order());
    }
    int order_themp = 0;
    for(_IT_GET_FEDATA_BY_TYPE_DOFS_FOR_LOOP_(this,termal_field_name,MBVERTEX,dit)) {
      order_themp = 1;
      break;
    }
    for(_IT_GET_FEDATA_BY_TYPE_DOFS_FOR_LOOP_(this,termal_field_name,MBEDGE,dit)) {
      order_themp = max(order_themp,dit->get_max_order());
    }
    int rule = (order+1)/2 + order_themp > 0 ? (order+1)/2 + order_themp : 0; //FIXME rule with temperature
    int nb_gauss_pts = gm_rule_size(rule,3);
    gaussPts.resize(4,nb_gauss_pts);
    ierr = Grundmann_Moeller_integration_points_3D_TET(
      rule,&gaussPts(0,0),&gaussPts(1,0),&gaussPts(2,0),&gaussPts(3,0)); CHKERRQ(ierr);
    g_NTET.resize(4*nb_gauss_pts);
    ierr = ShapeMBTET(&g_NTET[0],&gaussPts(0,0),&gaussPts(1,0),&gaussPts(2,0),nb_gauss_pts); CHKERRQ(ierr);
    g_TET_W = &(gaussPts(3,0));
    ierr = OpComplexForLazyStart(); CHKERRQ(ierr);
    ierr = GetIndicesMaterial(); CHKERRQ(ierr);
    ierr = GetData(
      order_x_edges,order_x_faces,order_x_volume,
      dofs_x_edge_data,dofs_x_edge,
      dofs_x_face_data,dofs_x_face,
      dofs_x_volume,dofs_x,
      spatial_field_name); CHKERRQ(ierr);
    switch(snes_ctx) {
      case CTX_SNESNONE:
      case CTX_SNESSETFUNCTION: { 
        ierr = GetFint(); CHKERRQ(ierr);
	ierr = AssembleMaterialFint(snes_f); CHKERRQ(ierr);
      }
      break;
      case CTX_SNESSETJACOBIAN: {
	ierr = GetTangent(); CHKERRQ(ierr);
	ierr = AssembleMaterialTangent(*snes_B); CHKERRQ(ierr);
	if(isCoupledProblem) {
	  ierr = GetIndicesRow(RowGlobSpatial,spatial_field_name); CHKERRQ(ierr);
	  ierr = AssembleMaterialCoupledTangent(*snes_B); CHKERRQ(ierr);
	}
      }
      break;
      default:
	SETERRQ(PETSC_COMM_SELF,1,"not implemented");
    }
    PetscFunctionReturn(0);
  }

};

struct MeshSmoothingFEMethod: public EshelbyFEMethod {

  MeshSmoothingFEMethod(FieldInterface& _mField,int _verbose = 0):
    FEMethod_ComplexForLazy_Data(_mField,_verbose),EshelbyFEMethod(_mField,0,0,_verbose) {
      type_of_analysis = mesh_quality_analysis;

      g_NTET.resize(4*1);
      ShapeMBTET(&g_NTET[0],G_TET_X1,G_TET_Y1,G_TET_Z1,1);
      g_TET_W = G_TET_W1;

    }

  PetscErrorCode preProcess() {
    PetscFunctionBegin;
    PetscFunctionReturn(0);
  }

  PetscErrorCode postProcess() {
    PetscFunctionBegin;
    PetscFunctionReturn(0);
  }

  virtual PetscErrorCode AssembleMeshSmoothingTangent(Mat B) {
    PetscFunctionBegin;
    switch(snes_ctx) {
      case CTX_SNESSETJACOBIAN:
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

  virtual PetscErrorCode AssembleMeshSmoothingFint(Vec f) {
    PetscFunctionBegin;
    switch(snes_ctx) {
      case CTX_SNESSETFUNCTION: { 
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

  PetscErrorCode operator()() {
    PetscFunctionBegin;
    ierr = OpComplexForLazyStart(); CHKERRQ(ierr);
    ierr = GetIndicesMaterial(); CHKERRQ(ierr);
    switch(snes_ctx) {
      case CTX_SNESSETFUNCTION:  
	ierr = GetFint(); CHKERRQ(ierr);
	ierr = AssembleMeshSmoothingFint(snes_f); CHKERRQ(ierr);
	break;
      case CTX_SNESSETJACOBIAN:
	ierr = GetTangent(); CHKERRQ(ierr);
	ierr = AssembleMeshSmoothingTangent(*snes_B); CHKERRQ(ierr);
	break;
      default:
	SETERRQ(PETSC_COMM_SELF,1,"not implemented");
    }
    PetscFunctionReturn(0);
  }

};

struct ArcLengthElemFEMethod: public FieldInterface::FEMethod {

  Interface& moab;
  ArcLengthCtx* arc_ptr;
  Vec GhostDiag;

  ArcLengthElemFEMethod(Interface& _moab,ArcLengthCtx *_arc_ptr): FEMethod(),moab(_moab),arc_ptr(_arc_ptr) {
    PetscInt ghosts[1] = { 0 };
    ParallelComm* pcomm = ParallelComm::get_pcomm(&moab,MYPCOMM_INDEX);
    if(pcomm->rank() == 0) {
      VecCreateGhost(PETSC_COMM_WORLD,1,1,0,ghosts,&GhostDiag);
    } else {
      VecCreateGhost(PETSC_COMM_WORLD,0,1,1,ghosts,&GhostDiag);
    }
  }
  ~ArcLengthElemFEMethod() {
    VecDestroy(&GhostDiag);
  }

  PetscErrorCode ierr;

  PetscErrorCode preProcess() {
    PetscFunctionBegin;
    switch(snes_ctx) {
      case CTX_SNESSETFUNCTION: { 
	ierr = calulate_dx_and_dlambda(snes_x); CHKERRQ(ierr);
	ierr = calulate_db(); CHKERRQ(ierr);
      }
      break;
      case CTX_SNESSETJACOBIAN: {
      }
      break;
      default:
      break;
    }
    PetscFunctionReturn(0);
  }

  double calulate_lambda_int() {
    PetscFunctionBegin;
    return arc_ptr->alpha*arc_ptr->dx2 + pow(arc_ptr->dlambda,2)*pow(arc_ptr->beta,2)*arc_ptr->F_lambda2;
    PetscFunctionReturn(0);
  }

  PetscErrorCode calulate_db() {
    PetscFunctionBegin;
    //db
    ierr = VecCopy(arc_ptr->dx,arc_ptr->db); CHKERRQ(ierr);
    ierr = VecScale(arc_ptr->db,2*arc_ptr->alpha); CHKERRQ(ierr);
    PetscFunctionReturn(0);
  }

  PetscErrorCode operator()() {
    PetscFunctionBegin;

    switch(snes_ctx) {
      case CTX_SNESSETFUNCTION: {
	arc_ptr->res_lambda = calulate_lambda_int() - pow(arc_ptr->s,2);
	ierr = VecSetValue(snes_f,arc_ptr->get_petsc_gloabl_dof_idx(),arc_ptr->res_lambda,ADD_VALUES); CHKERRQ(ierr);
	PetscPrintf(PETSC_COMM_SELF,"\tres_lambda = %6.4e\n",arc_ptr->res_lambda);
      }
      break; 
      case CTX_SNESSETJACOBIAN: {
	double diag = 2*arc_ptr->dlambda*pow(arc_ptr->beta,2)*arc_ptr->F_lambda2;
	ierr = VecSetValue(GhostDiag,0,diag,INSERT_VALUES); CHKERRQ(ierr);
	ierr = MatSetValue(*snes_B,arc_ptr->get_petsc_gloabl_dof_idx(),arc_ptr->get_petsc_gloabl_dof_idx(),1,ADD_VALUES); CHKERRQ(ierr);
      }
      break;
      default:
      break;
    }	
    
    PetscFunctionReturn(0);
  }

  PetscErrorCode postProcess() {
    PetscFunctionBegin;
    switch(snes_ctx) {
      case CTX_SNESSETFUNCTION: { 
	PetscPrintf(PETSC_COMM_WORLD,"\tlambda = %6.4e\n",arc_ptr->get_FieldData());  
      }
      break;
      case CTX_SNESSETJACOBIAN: {
	ierr = MatAssemblyBegin(*snes_B,MAT_FLUSH_ASSEMBLY); CHKERRQ(ierr);
	ierr = MatAssemblyEnd(*snes_B,MAT_FLUSH_ASSEMBLY); CHKERRQ(ierr);
	ierr = VecAssemblyBegin(GhostDiag); CHKERRQ(ierr);
	ierr = VecAssemblyEnd(GhostDiag); CHKERRQ(ierr);
	ierr = VecGhostUpdateBegin(GhostDiag,INSERT_VALUES,SCATTER_FORWARD); CHKERRQ(ierr);
	ierr = VecGhostUpdateEnd(GhostDiag,INSERT_VALUES,SCATTER_FORWARD); CHKERRQ(ierr);
	double *diag;
	ierr = VecGetArray(GhostDiag,&diag); CHKERRQ(ierr);
	arc_ptr->diag = *diag;
	ierr = VecRestoreArray(GhostDiag,&diag); CHKERRQ(ierr);
	PetscPrintf(PETSC_COMM_WORLD,"\tdiag = %6.4e\n",arc_ptr->diag);
      }
      break;
      default:
      break;
    }
    PetscFunctionReturn(0);
  }

  PetscErrorCode calulate_dx_and_dlambda(Vec x) {
    PetscFunctionBegin;
    //dx
    ierr = VecCopy(x,arc_ptr->dx); CHKERRQ(ierr);
    ierr = VecAXPY(arc_ptr->dx,-1,arc_ptr->x0); CHKERRQ(ierr);
    //dlambda
    if(arc_ptr->get_petsc_local_dof_idx()!=-1) {
      double *array;
      ierr = VecGetArray(arc_ptr->dx,&array); CHKERRQ(ierr);
      arc_ptr->dlambda = array[arc_ptr->get_petsc_local_dof_idx()];
      array[arc_ptr->get_petsc_local_dof_idx()] = 0;
      ierr = VecRestoreArray(arc_ptr->dx,&array); CHKERRQ(ierr);
    }
    int part = arc_ptr->get_part();
    MPI_Bcast(&(arc_ptr->dlambda),1,MPI_DOUBLE,part,PETSC_COMM_WORLD);
    //dx2
    ierr = VecDot(arc_ptr->dx,arc_ptr->dx,&arc_ptr->dx2); CHKERRQ(ierr);
    PetscPrintf(PETSC_COMM_WORLD,"\tdlambda = %6.4e dx2 = %6.4e\n",arc_ptr->dlambda,arc_ptr->dx2);
    PetscFunctionReturn(0);
  }

  PetscErrorCode calculate_init_dlambda(double *dlambda) {
    PetscFunctionBegin;
    *dlambda = sqrt(pow(arc_ptr->s,2)/(pow(arc_ptr->beta,2)*arc_ptr->F_lambda2));
    if(!(*dlambda == *dlambda)) {
      ostringstream sss;
      sss << "s " << arc_ptr->s << " " << arc_ptr->beta << " " << arc_ptr->F_lambda2;
      SETERRQ(PETSC_COMM_SELF,1,sss.str().c_str());
    }
    PetscFunctionReturn(0);
  }

  PetscErrorCode set_dlambda_to_x(Vec x,double dlambda) {
    PetscFunctionBegin;
    //check if locl dof idx is non zero, i.e. that lambda is acessible from this processor
    if(arc_ptr->get_petsc_local_dof_idx()!=-1) {
      double *array;
      ierr = VecGetArray(x,&array); CHKERRQ(ierr);
      double lambda_old = array[arc_ptr->get_petsc_local_dof_idx()];
      if(!(dlambda == dlambda)) {
	ostringstream sss;
	sss << "s " << arc_ptr->s << " " << arc_ptr->beta << " " << arc_ptr->F_lambda2;
	SETERRQ(PETSC_COMM_SELF,1,sss.str().c_str());
      }
      array[arc_ptr->get_petsc_local_dof_idx()] = lambda_old + dlambda;
      PetscPrintf(PETSC_COMM_WORLD,"\tlambda = %6.4e, %6.4e (%6.4e)\n",
	lambda_old, array[arc_ptr->get_petsc_local_dof_idx()], dlambda);
      ierr = VecRestoreArray(x,&array); CHKERRQ(ierr);
    }
    PetscFunctionReturn(0);
  }

};


}

#endif

