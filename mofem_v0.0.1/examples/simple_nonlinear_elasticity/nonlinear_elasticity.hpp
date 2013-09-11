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
#include "ElasticFEMethod.hpp"

#include "complex_for_lazy.h"

namespace MoFEM {

struct NL_ElasticFEMethod: public FEMethod_DriverComplexForLazy_Spatial {

  Range& SideSet1;
  Range& SideSet2;
  Range SideSet1_;

  NL_ElasticFEMethod(Interface& _moab,BaseDirihletBC *_dirihlet_bc_method_ptr,double _lambda,double _mu,Range &_SideSet1,Range &_SideSet2,int _verbose = 0): 
      FEMethod_DriverComplexForLazy_Spatial(_moab,_dirihlet_bc_method_ptr,_lambda,_mu,_verbose), SideSet1(_SideSet1),SideSet2(_SideSet2)  {

    set_PhysicalEquationNumber(neohookean);

    Range SideSet1Edges,SideSet1Nodes;
    rval = moab.get_adjacencies(SideSet1,1,false,SideSet1Edges,Interface::UNION); CHKERR_THROW(rval);
    rval = moab.get_connectivity(SideSet1,SideSet1Nodes,true); CHKERR_THROW(rval);
    SideSet1_.insert(SideSet1.begin(),SideSet1.end());
    SideSet1_.insert(SideSet1Edges.begin(),SideSet1Edges.end());
    SideSet1_.insert(SideSet1Nodes.begin(),SideSet1Nodes.end());

  }

  PetscErrorCode operator()() {
    PetscFunctionBegin;
    ierr = FEMethod_DriverComplexForLazy_Spatial::operator()(SideSet2); CHKERRQ(ierr);
    PetscFunctionReturn(0);
  }


};

}

#endif //__NONLINEAR_ELASTICITY_HPP__
