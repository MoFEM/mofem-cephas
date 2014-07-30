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

#ifndef __ELASTICFE_RVELAGRANGE_HOMOGENIZED_STRESS_PERIODIC_HPP__
#define __ELASTICFE_RVELAGRANGE_HOMOGENIZED_STRESS_PERIODIC_HPP__

#include "ElasticFE_RVELagrange_Periodic.hpp"

namespace MoFEM {
  
  struct ElasticFE_RVELagrange_Homogenized_Stress_Periodic: public ElasticFE_RVELagrange_Periodic {
    
    Vec DVec;
    Vec Stress_Homo;
    double *RVE_volume;
    
    ElasticFE_RVELagrange_Homogenized_Stress_Periodic(FieldInterface& _mField,Mat &_Aij,Vec &_D,Vec& _F,double *_RVE_volume, ublas::vector<FieldData> _applied_strain, Vec& _Stress_Homo):
    ElasticFE_RVELagrange_Periodic(_mField,_Aij, _D, _F, _applied_strain), DVec(_D), Stress_Homo(_Stress_Homo),RVE_volume(_RVE_volume) {};
    
    PetscErrorCode preProcess();
    PetscErrorCode postProcess();
    
    ublas::vector<ublas::vector<FieldData> > Lamda;
    virtual PetscErrorCode Calculate_Homo_Stress();
    PetscErrorCode operator()();
  };
  
}

#endif //__ElasticFE_RVELagrange_Periodic_RigidBodyMotion
