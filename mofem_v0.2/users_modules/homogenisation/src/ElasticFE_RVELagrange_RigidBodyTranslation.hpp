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

#ifndef __ELASTICFE_RVELAGRANGE_RIGIDBODYTRANSLATION_HPP__
#define __ELASTICFE_RVELAGRANGE_RIGIDBODYTRANSLATION_HPP__


using namespace ObosleteUsersModules;

namespace MoFEM {
  
  struct ElasticFE_RVELagrange_RigidBodyTranslation: public ElasticFE_RVELagrange_Disp {
    
    const string field_lagrange_rigid_tans;
    
    ElasticFE_RVELagrange_RigidBodyTranslation(FieldInterface& _mField,Mat &_Aij,Vec &_D,Vec& _F, ublas::vector<FieldData> _applied_strain, const string& _field_main, const string& _field_lagrange, int _rank_field, const string& _field_lagrange_rigid_tans):
    ElasticFE_RVELagrange_Disp(_mField,_Aij, _D, _F, _applied_strain, _field_main, _field_lagrange, _rank_field), field_lagrange_rigid_tans(_field_lagrange_rigid_tans){};
    
    
    virtual PetscErrorCode GetN_and_Indices();
    virtual PetscErrorCode Lhs();
    PetscErrorCode operator()();
    
  };
}

#endif //__ElasticFE_RVELagrange_Periodic_RigidBodyMotion
