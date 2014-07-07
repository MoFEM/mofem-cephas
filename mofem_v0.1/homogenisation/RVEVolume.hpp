/* Copyright (C) 2014, Zahur Ullah (Zahur.Ullah AT glasgow.ac.uk)
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

#ifndef __RVEVolume_HPP__
#define __RVEVolume_HPP__

#include <boost/numeric/ublas/symmetric.hpp>
#include "ElasticFEMethod.hpp"
#include "ElasticFEMethodTransIso.hpp"


namespace MoFEM {
  
  struct RVEVolume: public ElasticFEMethod {
    Vec RVE_volume_Vec;
    
    RVEVolume(FieldInterface& _mField,BaseDirihletBC *_dirihlet_ptr,Mat &_Aij,Vec &_D,Vec& _F,double _lambda,double _mu, Vec _RVE_volume_Vec):
    ElasticFEMethod(_mField,_dirihlet_ptr,_Aij,_D,_F,_lambda,_mu), RVE_volume_Vec(_RVE_volume_Vec){};
    
    
    PetscErrorCode postProcess() {
      PetscFunctionBegin;
      ierr = VecAssemblyBegin(RVE_volume_Vec); CHKERRQ(ierr);
      ierr = VecAssemblyEnd(RVE_volume_Vec); CHKERRQ(ierr);
      PetscFunctionReturn(0);
    }
    
    PetscErrorCode operator()() {
      PetscFunctionBegin;
      ierr = OpStudentStart_TET(g_NTET); CHKERRQ(ierr);
      
      int Indices[1];  Indices[0]=pcomm->rank();
      double Vol_elm[1];  Vol_elm[0]=V;
      ierr = VecSetValues(RVE_volume_Vec,1,Indices,Vol_elm,ADD_VALUES); CHKERRQ(ierr);
      
      ierr = OpStudentEnd(); CHKERRQ(ierr);
      PetscFunctionReturn(0);
    }
  };
  
  
  
  struct RVEVolumeTrans: public TranIsotropicFibreDirRotElasticFEMethod {
    Vec RVE_volume_Vec;
    
    RVEVolumeTrans(FieldInterface& _mField,BaseDirihletBC *_dirihlet_ptr,Mat &_Aij,Vec &_D,Vec& _F, Vec _RVE_volume_Vec):
    TranIsotropicFibreDirRotElasticFEMethod(_mField,_dirihlet_ptr,_Aij,_D,_F), RVE_volume_Vec(_RVE_volume_Vec){};
    
    
    PetscErrorCode postProcess() {
      PetscFunctionBegin;
      ierr = VecAssemblyBegin(RVE_volume_Vec); CHKERRQ(ierr);
      ierr = VecAssemblyEnd(RVE_volume_Vec); CHKERRQ(ierr);
      PetscFunctionReturn(0);
    }
    
    PetscErrorCode operator()() {
      PetscFunctionBegin;
      ierr = OpStudentStart_TET(g_NTET); CHKERRQ(ierr);
      
      int Indices[1];  Indices[0]=pcomm->rank();
      double Vol_elm[1];  Vol_elm[0]=V;
      ierr = VecSetValues(RVE_volume_Vec,1,Indices,Vol_elm,ADD_VALUES); CHKERRQ(ierr);
      
      ierr = OpStudentEnd(); CHKERRQ(ierr);
      PetscFunctionReturn(0);
    }
  };
  
  
  
  
  
  
}

#endif //__RVEVolume_HPP__
