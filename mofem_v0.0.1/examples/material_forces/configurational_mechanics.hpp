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

#include "moabField.hpp"

using namespace MoFEM;

struct ConfigurationalMechanics {
 
  Tag th_MaterialFireWall;
  typedef bitset<16> Material_FirelWall_def;
  Material_FirelWall_def *Material_FirelWall;

  EntityHandle CornersNodesMeshset,SurfacesFacesMeshset,CrackSurfacesFacesMeshset,CrackForntMeshset;
  
  PetscErrorCode ConfigurationalMechanics_SetMaterialFireWall(moabField& mField);
  PetscErrorCode ConfigurationalMechanics_SpatialProblemDefinition(moabField& mField); 
  PetscErrorCode ConfigurationalMechanics_MaterialProblemDefinition(moabField& mField);
  PetscErrorCode ConfigurationalMechanics_CoupledProblemDefinition(moabField& mField);
  PetscErrorCode ConfigurationalMechanics_ConstrainsProblemDefinition(moabField& mField);
  PetscErrorCode ConfigurationalMechanics_ConstrainsCrackFrontProblemDefinition(moabField& mField);
  PetscErrorCode ConfigurationalMechanics_SpatialPartitionProblems(moabField& mField);
  PetscErrorCode ConfigurationalMechanics_MaterialPartitionProblems(moabField& mField);
  PetscErrorCode ConfigurationalMechanics_CoupledPartitionProblems(moabField& mField);
  PetscErrorCode ConfigurationalMechanics_ConstrainsPartitionProblems(moabField& mField,string problem);
  PetscErrorCode ConfigurationalMechanics_CrackFrontPartitionProblems(moabField& mField,string problem);
  PetscErrorCode ConfigurationalMechanics_SetSpatialPositions(moabField& mField);
  PetscErrorCode ConfigurationalMechanics_SetMaterialPositions(moabField& mField);
  PetscErrorCode ConfigurationalMechanics_SolveSpatialProblem(moabField& mField);
  PetscErrorCode ConfigurationalMechanics_SolveMaterialProblem(moabField& mField);
  PetscErrorCode ConfigurationalMechanics_CalculateSpatialResidual(moabField& mField);
  PetscErrorCode ConfigurationalMechanics_CalculateMaterialForces(moabField& mField,string problem);
  PetscErrorCode ConfigurationalMechanics_ProjectForceVector(moabField& mField,string problem);
  PetscErrorCode ConfigurationalMechanics_GriffithForceVector(moabField& mField);
  PetscErrorCode ConfigurationalMechanics_SolveCoupledProblem(moabField& mField);

};

#endif //__CONFIGURATIONAL_MECHANICS_HPP__
