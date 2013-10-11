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

#ifndef __CONFIGURATIONAL_MECHANICS_HPP__
#define __CONFIGURATIONAL_MECHANICS_HPP__

#include "FieldInterface.hpp"
#include "petscShellMATs_ConstrainsByMarkAinsworth.hpp"

using namespace MoFEM;

struct ConfigurationalMechanics {
 
  Tag th_MaterialFireWall;
  typedef bitset<16> Material_FirelWall_def;
  Material_FirelWall_def *material_FirelWall;

  enum FirWall {
    FW_spatial_problem_definition = 0,
    FW_material_problem_definition,
    FW_coupled_problem_definition,
    FW_constrains_problem_definition,
    FW_constrains_crack_front_problem_definition,
    FW_set_spatial_positions,
    FW_set_material_positions
  };

  EntityHandle cornersNodesMeshset,surfacesFacesMeshset,cracksurfacesFacesMeshset,crackForntMeshset;
  matPROJ_ctx *projSurfaceCtx,*projFrontCtx;
  
  PetscErrorCode set_material_fire_wall(FieldInterface& mField);
  PetscErrorCode spatial_problem_definition(FieldInterface& mField); 
  PetscErrorCode material_problem_definition(FieldInterface& mField);
  PetscErrorCode coupled_problem_definition(FieldInterface& mField);
  PetscErrorCode constrains_problem_definition(FieldInterface& mField);
  PetscErrorCode constrains_crack_front_problem_definition(FieldInterface& mField);
  PetscErrorCode spatialPartitionProblems(FieldInterface& mField);
  PetscErrorCode material_partition_problems(FieldInterface& mField);
  PetscErrorCode coupled_partition_problems(FieldInterface& mField);
  PetscErrorCode constrains_partition_problems(FieldInterface& mField,string problem);
  PetscErrorCode crackfront_partition_problems(FieldInterface& mField,string problem);
  PetscErrorCode set_spatial_positions(FieldInterface& mField);
  PetscErrorCode set_material_positions(FieldInterface& mField);
  PetscErrorCode solve_spatial_problem(FieldInterface& mField);
  PetscErrorCode solve_material_problem(FieldInterface& mField);
  PetscErrorCode calculate_spatial_residual(FieldInterface& mField);
  PetscErrorCode calculate_material_forces(FieldInterface& mField,string problem);
  PetscErrorCode surface_projection_data(FieldInterface& mField,string problem);
  PetscErrorCode project_force_vector(FieldInterface& mField,string problem);
  PetscErrorCode front_projection_data(FieldInterface& mField,string problem);
  PetscErrorCode griffith_force_vector(FieldInterface& mField,string problem);
  PetscErrorCode griffith_g(FieldInterface& mField,string problem);
  PetscErrorCode solve_coupled_problem(FieldInterface& mField);

};

#endif //__CONFIGURATIONAL_MECHANICS_HPP__
