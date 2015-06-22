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

#ifndef __MAIN_CRACK_FUNCTION__
#define __MAIN_CRACK_FUNCTION__

PetscErrorCode main_spatial_solution(FieldInterface& m_field,ConfigurationalFractureMechanics& conf_prob);
PetscErrorCode main_material_forces(FieldInterface& m_field,ConfigurationalFractureMechanics& conf_prob);

//crack propagation

/** \brief rescale load factor, such that maximally stressed crack front node has griffithe energy equal to gc
  *
  */
PetscErrorCode main_rescale_load_factor(FieldInterface& m_field,ConfigurationalFractureMechanics& conf_prob);

PetscErrorCode main_arc_length_setup(FieldInterface& m_field,ConfigurationalFractureMechanics& conf_prob);
PetscErrorCode main_arc_length_restart(FieldInterface& m_field,ConfigurationalFractureMechanics& conf_prob);
PetscErrorCode main_arc_length_solve(FieldInterface& m_field,ConfigurationalFractureMechanics& conf_prob);

#endif //__MAIN_CRACK_FUNCTION__
