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

#ifndef __NONLINEAR_ELASTICITY_HPP__
#define __NONLINEAR_ELASTICITY_HPP__

#include "moabField.hpp"
#include "moabField_Core.hpp"
#include "moabFEMethod_UpLevelStudent.hpp"
#include "cholesky.hpp"
#include <petscksp.h>

#include "moabSnes.hpp"
#include "moabFEMethod_ComplexForLazy.hpp"
#include "moabFEMethod_DriverComplexForLazy.hpp"
#include "moabConstrainsByMarkAinsworth.hpp"
#include "moabSurfaceConstrains.hpp"
#include "ElasticFEMethod.hpp"

#include "complex_for_lazy.h"

namespace MoFEM {

struct SetPositionsEntMethod: public moabField::EntMethod {
    ErrorCode rval;
    PetscErrorCode ierr;
    Interface& moab;

    EntityHandle node;
    double coords[3];

    SetPositionsEntMethod(Interface& _moab): EntMethod(),moab(_moab),node(no_handle) {}
    
    PetscErrorCode preProcess() {
      PetscFunctionBegin;
      PetscPrintf(PETSC_COMM_WORLD,"Start Set Positions\n");
      PetscFunctionReturn(0);
    } 
     
    PetscErrorCode operator()() {
      PetscFunctionBegin;
      if(dof_ptr->get_ent_type()!=MBVERTEX) PetscFunctionReturn(0);
      EntityHandle ent = dof_ptr->get_ent();
      int dof_rank = dof_ptr->get_dof_rank();
      double &fval = dof_ptr->get_FieldData();
      if(node!=ent) {
	rval = moab.get_coords(&ent,1,coords); CHKERR_PETSC(rval);
	node = ent;
      }
      fval = coords[dof_rank];
      PetscFunctionReturn(0);
    }

    PetscErrorCode postProcess() {
      PetscFunctionBegin;
      PetscPrintf(PETSC_COMM_WORLD,"End Set Positions\n");
      PetscFunctionReturn(0);
    }

};

struct Spatial_ElasticFEMethod: public FEMethod_DriverComplexForLazy_Spatial {

  Range& SideSet2;

  Spatial_ElasticFEMethod(Interface& _moab,BaseDirihletBC *_dirihlet_bc_method_ptr,double _lambda,double _mu,Range &_SideSet2,int _verbose = 0): 
      FEMethod_DriverComplexForLazy_Spatial(_moab,_dirihlet_bc_method_ptr,_lambda,_mu,_verbose), SideSet2(_SideSet2)  {

    set_PhysicalEquationNumber(hooke);
    //set_PhysicalEquationNumber(neohookean);

  }

  PetscErrorCode operator()() {
    PetscFunctionBegin;
    ierr = FEMethod_DriverComplexForLazy_Spatial::operator()(SideSet2); CHKERRQ(ierr);
    PetscFunctionReturn(0);
  }


};

struct MaterialForcesFEMethod: public FEMethod_DriverComplexForLazy_Material {

  Vec F_MATERIAL;
  MaterialForcesFEMethod(Interface& _moab,BaseDirihletBC *_dirihlet_bc_method_ptr,double _lambda,double _mu,Vec _F_MATERIAL,int _verbose = 0): 
      FEMethod_DriverComplexForLazy_Material(_moab,_dirihlet_bc_method_ptr,_lambda,_mu,_verbose),F_MATERIAL(_F_MATERIAL) {

    set_PhysicalEquationNumber(hooke);
    type_of_analysis = material_analysis;

  }

  PetscErrorCode preProcess() {
    PetscFunctionBegin;
    PetscFunctionReturn(0);
  }

  PetscErrorCode operator()() {
    PetscFunctionBegin;
    ierr = OpComplexForLazyStart(); CHKERRQ(ierr);
    ierr = GetIndices(RowGlobMaterial,ColGlobMaterial,material_field_name); CHKERRQ(ierr);
    ierr = GetData(dofs_x_edge_data,dofs_x_edge,
      dofs_x_face_data,dofs_x_face,
      dofs_x_volume,dofs_x,
      spatial_field_name); CHKERRQ(ierr);
    ierr = GetFint(); CHKERRQ(ierr);
    ierr = VecSetValues(F_MATERIAL,RowGlobMaterial[0].size(),&(RowGlobMaterial[0])[0],&(Fint_H.data())[0],ADD_VALUES); CHKERRQ(ierr);
    PetscFunctionReturn(0);
  }

  PetscErrorCode postProcess() {
    PetscFunctionBegin;
    PetscFunctionReturn(0);
  }

};

struct Material_ElasticFEMethod: public FEMethod_DriverComplexForLazy_Material {

  moabField& mField;
  Material_ElasticFEMethod(Interface& _moab,moabField& _mField,BaseDirihletBC *_dirihlet_bc_method_ptr,double _lambda,double _mu,int _verbose = 0): 
      FEMethod_DriverComplexForLazy_Material(_moab,_dirihlet_bc_method_ptr,_lambda,_mu,_verbose),mField(_mField)  {
    set_PhysicalEquationNumber(hooke);
    //set_PhysicalEquationNumber(neohookean);
  }

};


}

#endif //__NONLINEAR_ELASTICITY_HPP__
