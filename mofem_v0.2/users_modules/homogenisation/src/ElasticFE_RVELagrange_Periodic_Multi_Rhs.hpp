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

#ifndef __ELASTICFE_RVELAGRANGE_PERIODIC_MULTI_RHS_HPP__
#define __ELASTICFE_RVELAGRANGE_PERIODIC_MULTI_RHS_HPP__

using namespace ObosleteUsersModules;


namespace MoFEM {
  
  
  struct ElasticFE_RVELagrange_Periodic_Multi_Rhs: public ElasticFE_RVELagrange_Periodic {
    
    ublas::vector<FieldData> _applied_strain;
    
    Vec F1,F2,F3,F4,F5,F6;
    
    
    //constructor for mechanical problem (rank=3)
    ElasticFE_RVELagrange_Periodic_Multi_Rhs(
                                             FieldInterface& _mField,Mat &_Aij,Vec &_D,Vec& _F1,Vec& _F2,Vec& _F3,Vec& _F4,Vec& _F5,Vec& _F6,ublas::vector<FieldData> _applied_strain, const string& _field_main, const string& _field_lagrange, int _rank_field):
    ElasticFE_RVELagrange_Periodic(_mField,_Aij, _D, _F1, _applied_strain, _field_main, _field_lagrange, _rank_field),F1(_F1),F2(_F2),F3(_F3),F4(_F4),F5(_F5),F6(_F6){};
    
    //constructor for thermal and moisture transport problems (rank=1)
    ElasticFE_RVELagrange_Periodic_Multi_Rhs(
                                             FieldInterface& _mField,Mat &_Aij,Vec &_D,Vec& _F1,Vec& _F2,Vec& _F3,ublas::vector<FieldData> _applied_strain, const string& _field_main, const string& _field_lagrange, int _rank_field):
    ElasticFE_RVELagrange_Periodic(_mField,_Aij, _D, _F1, _applied_strain, _field_main, _field_lagrange, _rank_field),F1(_F1),F2(_F2),F3(_F3){};
    
    virtual PetscErrorCode postProcess();
    virtual PetscErrorCode Rhs();
    
    
  };

  
}

#endif //__ElasticFE_RVELagrange_Periodic__
