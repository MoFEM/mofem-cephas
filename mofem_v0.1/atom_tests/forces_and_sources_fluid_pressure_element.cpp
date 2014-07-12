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
#include "FluidPressure.hpp"

#include <boost/iostreams/tee.hpp>
#include <boost/iostreams/stream.hpp>
#include <fstream>
#include <iostream>

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

  PetscBool flg = PETSC_TRUE;
  char mesh_file_name[255];
  ierr = PetscOptionsGetString(PETSC_NULL,"-my_file",mesh_file_name,255,&flg); CHKERRQ(ierr);
  if(flg != PETSC_TRUE) {
    SETERRQ(PETSC_COMM_SELF,1,"*** ERROR -my_file (MESH FILE NEEDED)");
  }

  const char *option;
  option = "";//"PARALLEL=BCAST;";//;DEBUG_IO";
  rval = moab.load_file(mesh_file_name, 0, option); CHKERR_PETSC(rval); 

  //Create MoFEM (Joseph) database
  FieldCore core(moab);
  FieldInterface& mField = core;

  //set entitities bit level
  BitRefLevel bit_level0;
  bit_level0.set(0);
  EntityHandle meshset_level0;
  rval = moab.create_meshset(MESHSET_SET,meshset_level0); CHKERR_PETSC(rval);
  ierr = mField.seed_ref_level_3D(0,bit_level0); CHKERRQ(ierr);

  //Definitions

  //add DISPLACEMENT field, Hilbert space H1, veror field rank 3 (displacemnt
  //has three components ux,uy,uz)
  ierr = mField.add_field("DISPLACEMENT",H1,3); CHKERRQ(ierr);

  ///add probelem which will be solved, could be more than one problem
  //operating on some subset of defined approximatons spces
  ierr = mField.add_problem("TEST_PROBLEM"); CHKERRQ(ierr);

  //mesh could have several refinment levels which share some subset of entities between them. 
  //below defines on which set of entities (on refinment level 0) build approximation spaces for TEST_PROBLEM
  ierr = mField.modify_problem_ref_level_add_bit("TEST_PROBLEM",bit_level0); CHKERRQ(ierr);

  //add entities on wich DISPLACEMENT field is approximated, you can add
  //entities form several approximation levels at onec. You can as well
  //approximate field only on some mesh subdomain, in that case displacements
  //are approximated on roor moab mesh.
  ierr = mField.add_ents_to_field_by_TETs(0,"DISPLACEMENT"); CHKERRQ(ierr);

  //set app. order for displacement fiedl. it is set uniform approximation
  //order. in genreal evry entity can have arbitraty approximation level,
  //ranging from 1 to 10 and more.
  int order = 1;
  ierr = mField.set_field_order(0,MBTET,"DISPLACEMENT",order); CHKERRQ(ierr);
  ierr = mField.set_field_order(0,MBTRI,"DISPLACEMENT",order); CHKERRQ(ierr);
  ierr = mField.set_field_order(0,MBEDGE,"DISPLACEMENT",order); CHKERRQ(ierr);
  ierr = mField.set_field_order(0,MBVERTEX,"DISPLACEMENT",1); CHKERRQ(ierr);

  //define fluid pressure finite elements
  FluidPressure fluid_pressure_fe(mField);
  fluid_pressure_fe.addNeumannFluidPressureBCElements("TEST_PROBLEM","DISPLACEMENT");

  //construct data structrures for fields and finite elements. at that points
  //entities, finite elements or dofs have unque uid, but are not partitioned
  //or numbered. user can add entities to mesh, add dofs or elenents if
  //necessaery. in case of modifications data structures are updated.
  ierr = mField.build_fields(); CHKERRQ(ierr);
  ierr = mField.build_finiteElementsPtr(); CHKERRQ(ierr);
  ierr = mField.build_adjacencies(bit_level0); CHKERRQ(ierr);
  ierr = mField.build_problems(); CHKERRQ(ierr);

  //to solve problem it needt to be respresented in matrix vector form. this
  //demand numbertion of dofs and proble  partitioning.
  ierr = mField.simple_partition_problem("TEST_PROBLEM"); CHKERRQ(ierr);
  ierr = mField.partition_finiteElementsPtr("TEST_PROBLEM"); CHKERRQ(ierr);
  //what are ghost nodes, see Petsc Manual
  ierr = mField.partition_ghost_dofs("TEST_PROBLEM"); CHKERRQ(ierr);

  //create vector for problem
  Vec F;
  ierr = mField.VecCreateGhost("TEST_PROBLEM",ROW,&F); CHKERRQ(ierr);
  ierr = fluid_pressure_fe.setNeumannFluidPressureFiniteElementOperators("DISPLACEMENT",F,false,false); CHKERRQ(ierr);

  ierr = VecZeroEntries(F); CHKERRQ(ierr);
  ierr = mField.loop_finite_elements("TEST_PROBLEM","FLUID_PRESSURE_FE",fluid_pressure_fe.getLoopFe()); CHKERRQ(ierr);
  ierr = VecAssemblyBegin(F); CHKERRQ(ierr);
  ierr = VecAssemblyEnd(F); CHKERRQ(ierr);
  ierr = mField.set_global_VecCreateGhost("TEST_PROBLEM",ROW,F,INSERT_VALUES,SCATTER_REVERSE); CHKERRQ(ierr);

  //ierr = VecView(F,PETSC_VIEWER_STDOUT_WORLD); CHKERRQ(ierr);

  PetscViewer viewer;
  ierr = PetscViewerASCIIOpen(PETSC_COMM_WORLD,"forces_and_sources_fluid_pressure_element.txt",&viewer); CHKERRQ(ierr);
  ierr = VecChop(F,1e-4); CHKERRQ(ierr);
  ierr = VecView(F,viewer); CHKERRQ(ierr);
  ierr = PetscViewerDestroy(&viewer); CHKERRQ(ierr);

  /*double sum = 0;
  ierr = VecSum(F,&sum); CHKERRQ(ierr);
  ierr = PetscPrintf(PETSC_COMM_WORLD,"sum  = %4.3f\n",sum); CHKERRQ(ierr);

  map<EntityHandle,ublas::vector<double> > tags_vals;
  for(_IT_GET_DOFS_FIELD_BY_NAME_FOR_LOOP_(mField,"DISPLACEMENT",dof)) {
    tags_vals[dof->get_ent()].resize(3);
    tags_vals[dof->get_ent()][dof->get_dof_rank()] = dof->get_FieldData();
  }
  vector<EntityHandle> ents;
  ents.resize(tags_vals.size());
  vector<double> vals(3*tags_vals.size());
  int idx = 0;
  for(map<EntityHandle,ublas::vector<double> >::iterator mit = tags_vals.begin();
    mit!=tags_vals.end();mit++,idx++) {
    ents[idx] = mit->first;
    vals[3*idx + 0] = mit->second[0];
    vals[3*idx + 1] = mit->second[1];   
    vals[3*idx + 2] = mit->second[2];
  }
  
  double def_VAL[3] = {0,0,0};
  Tag th_vals;
  rval = moab.tag_get_handle("FLUID_PRESURE_FORCES",3,MB_TYPE_DOUBLE,th_vals,MB_TAG_CREAT|MB_TAG_SPARSE,def_VAL); CHKERR_PETSC(rval);
  rval = moab.tag_set_data(th_vals,&ents[0],ents.size(),&vals[0]); CHKERR_PETSC(rval);

  EntityHandle out_meshset;
  rval = moab.create_meshset(MESHSET_SET,out_meshset); CHKERR_PETSC(rval);
  ierr = mField.problem_get_FE("TEST_PROBLEM","FLUID_PRESSURE_FE",out_meshset); CHKERRQ(ierr);
  rval = moab.write_file("out.vtk","VTK","",&out_meshset,1); CHKERR_PETSC(rval);
  rval = moab.delete_entities(&out_meshset,1); CHKERR_PETSC(rval);*/
  
  //destroy vector
  ierr = VecDestroy(&F); CHKERRQ(ierr);

  ierr = PetscFinalize(); CHKERRQ(ierr);

  return 0;

}


