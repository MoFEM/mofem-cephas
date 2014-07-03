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

#include "FieldInterface.hpp"
#include "FieldCore.hpp"
#include "ThermalStressElement.hpp"

#include "PostProcVertexMethod.hpp"
#include "PostProcDisplacementAndStrainOnRefindedMesh.hpp"

#include <boost/iostreams/tee.hpp>
#include <boost/iostreams/stream.hpp>
#include <fstream>
#include <iostream>

#include <petscksp.h>

namespace bio = boost::iostreams;
using bio::tee_device;
using bio::stream;

using namespace MoFEM;

static char help[] = "...\n\n";

int main(int argc, char *argv[]) {

  ErrorCode rval;
  PetscErrorCode ierr;

  PetscInitialize(&argc,&argv,(char *)0,help);

  Core mb_instance;
  Interface& moab = mb_instance;
  int rank;
  MPI_Comm_rank(PETSC_COMM_WORLD,&rank);

  PetscBool flg = PETSC_TRUE;
  char mesh_file_name[255];
  ierr = PetscOptionsGetString(PETSC_NULL,"-my_file",mesh_file_name,255,&flg); CHKERRQ(ierr);
  if(flg != PETSC_TRUE) {
    SETERRQ(PETSC_COMM_SELF,1,"*** ERROR -my_file (MESH FILE NEEDED)");
  }

  ParallelComm* pcomm = ParallelComm::get_pcomm(&moab,MYPCOMM_INDEX);
  if(pcomm == NULL) pcomm =  new ParallelComm(&moab,PETSC_COMM_WORLD);

  const char *option;
  option = "";//"PARALLEL=BCAST;";//;DEBUG_IO";
  BARRIER_RANK_START(pcomm) 
  rval = moab.load_file(mesh_file_name, 0, option); CHKERR_PETSC(rval); 
  BARRIER_RANK_END(pcomm) 

  //Create MoFEM (Joseph) database
  FieldCore core(moab);
  FieldInterface& mField = core;

  //set entitities bit level
  BitRefLevel bit_level0;
  bit_level0.set(0);
  EntityHandle meshset_level0;
  rval = moab.create_meshset(MESHSET_SET,meshset_level0); CHKERR_PETSC(rval);
  ierr = mField.seed_ref_level_3D(0,bit_level0); CHKERRQ(ierr);

  //Fields
  ierr = mField.add_field("DISP",H1,3); CHKERRQ(ierr);
  ierr = mField.add_field("TEMP",H1,1); CHKERRQ(ierr);

  //Problem
  ierr = mField.add_problem("PROB"); CHKERRQ(ierr);

  //set refinment level for problem
  ierr = mField.modify_problem_ref_level_add_bit("PROB",bit_level0); CHKERRQ(ierr);

  //meshset consisting all entities in mesh
  EntityHandle root_set = moab.get_root_set(); 
  //add entities to field
  ierr = mField.add_ents_to_field_by_TETs(root_set,"TEMP"); CHKERRQ(ierr);
  ierr = mField.add_ents_to_field_by_TETs(root_set,"DISP"); CHKERRQ(ierr);

  //set app. order
  //see Hierarchic Finite Element Bases on Unstructured Tetrahedral Meshes (Mark Ainsworth & Joe Coyle)
  int order = 4;
  ierr = mField.set_field_order(root_set,MBTET,"TEMP",order); CHKERRQ(ierr);
  ierr = mField.set_field_order(root_set,MBTRI,"TEMP",order); CHKERRQ(ierr);
  ierr = mField.set_field_order(root_set,MBEDGE,"TEMP",order); CHKERRQ(ierr);
  ierr = mField.set_field_order(root_set,MBVERTEX,"TEMP",1); CHKERRQ(ierr);

  ierr = mField.set_field_order(root_set,MBTET,"DISP",order); CHKERRQ(ierr);
  ierr = mField.set_field_order(root_set,MBTRI,"DISP",order); CHKERRQ(ierr);
  ierr = mField.set_field_order(root_set,MBEDGE,"DISP",order); CHKERRQ(ierr);
  ierr = mField.set_field_order(root_set,MBVERTEX,"DISP",1); CHKERRQ(ierr);

  ThermalStressElement thermal_stress_elem(mField);
  ierr = thermal_stress_elem.addThermalSterssElement("PROB","ELAS","DISP","TEMP"); CHKERRQ(ierr);

  /****/
  //build database
  //build field
  ierr = mField.build_fields(); CHKERRQ(ierr);
  //build finite elemnts
  ierr = mField.build_finite_elements(); CHKERRQ(ierr);
  //build adjacencies
  ierr = mField.build_adjacencies(bit_level0); CHKERRQ(ierr);
  //build problem
  ierr = mField.build_problems(); CHKERRQ(ierr);

  /****/
  //mesh partitioning 
  //partition
  ierr = mField.simple_partition_problem("TEST_PROBLEM"); CHKERRQ(ierr);
  ierr = mField.partition_finite_elements("TEST_PROBLEM"); CHKERRQ(ierr);
  //what are ghost nodes, see Petsc Manual
  ierr = mField.partition_ghost_dofs("TEST_PROBLEM"); CHKERRQ(ierr);

  

  ierr = PetscFinalize(); CHKERRQ(ierr);

  return 0;

}


