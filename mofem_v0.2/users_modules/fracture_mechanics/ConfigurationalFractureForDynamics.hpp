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

#ifndef __CONFIGURATIONAL_DYNAMICS_HPP__
#define __CONFIGURATIONAL_DYNAMICS_HPP__

using namespace ObosleteUsersModules;

struct ConfigurationalFracturDynamics: public ConfigurationalFractureMechanics {

  ConfigurationalFracturDynamics(FieldInterface& m_field):
    ConfigurationalFractureMechanics(m_field) {}

  PetscErrorCode coupled_dynamic_problem_definition(FieldInterface& m_field);
  PetscErrorCode coupled_dynamic_partition_problems(FieldInterface& m_field);
  PetscErrorCode solve_dynmaic_problem(FieldInterface& m_field,TS *ts,double fraction_treshold = 1e-1);

}; 

#endif //__CONFIGURATIONAL_DYNAMICS_HPP__


