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

#ifndef __ElasticFE_RVELagrange_RigidBodyTranslation_HPP__
#define __ElasticFE_RVELagrange_RigidBodyTranslation__

#include "ElasticFE_RVELagrange_Disp.hpp"

namespace MoFEM {
  
  struct ElasticFE_RVELagrange_RigidBodyTranslation: public ElasticFE_RVELagrange_Disp {
    
    ElasticFE_RVELagrange_RigidBodyTranslation(FieldInterface& _mField,Mat &_Aij,Vec &_D,Vec& _F, ublas::vector<FieldData> _applied_strain):
    ElasticFE_RVELagrange_Disp(_mField,_Aij, _D, _F, _applied_strain){};
    
    virtual PetscErrorCode GetN_and_Indices();
    virtual PetscErrorCode Lhs();
    PetscErrorCode operator()();
    
  };
}

#endif //__ElasticFE_RVELagrange_Periodic_RigidBodyMotion
