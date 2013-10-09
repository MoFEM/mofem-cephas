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

using namespace MoFEM;

ErrorCode rval;
PetscErrorCode ierr;

static char help[] = "teting mesh refinment algorithm\n\n";

int main(int argc, char *argv[]) {

  PetscInitialize(&argc,&argv,(char *)0,help);

  PetscBool flg = PETSC_TRUE;
  char mesh_file_name[255];
  ierr = PetscOptionsGetString(PETSC_NULL,"-my_file",mesh_file_name,255,&flg); CHKERRQ(ierr);
  if(flg != PETSC_TRUE) {
    SETERRQ(PETSC_COMM_SELF,1,"*** ERROR -my_file (MESH FILE NEEDED)");
  }
  char mesh_out_file_name[255];
  ierr = PetscOptionsGetString(PETSC_NULL,"-my_out_file",mesh_out_file_name,255,&flg); CHKERRQ(ierr);
  if(flg != PETSC_TRUE) {
    SETERRQ(PETSC_COMM_SELF,1,"*** ERROR -my_out_file (MESH FILE NEEDED)");
  }

  Core mb_instance;
  Interface& moab = mb_instance;

  const char *option;
  option = "";//"PARALLEL=BCAST";//;DEBUG_IO";
  rval = moab.load_file(mesh_file_name, 0, option); CHKERR(rval); 
  ParallelComm* pcomm = ParallelComm::get_pcomm(&moab,0);

  PetscLogDouble t1,t2;
  PetscLogDouble v1,v2;
  ierr = PetscGetTime(&v1); CHKERRQ(ierr);
  ierr = PetscGetCPUTime(&t1); CHKERRQ(ierr);

  FieldCore core(moab);
  FieldInterface& mField = core;

  Range CubitSideSets_meshsets;
  ierr = mField.get_CubitBCType_meshsets(SideSet,CubitSideSets_meshsets); CHKERRQ(ierr);

  ierr = mField.seed_ref_level_3D(0,0); CHKERRQ(ierr);
  BitRefLevel bit_level0;
  bit_level0.set(0);
  Range::iterator mit = CubitSideSets_meshsets.begin();
  for(;mit!=CubitSideSets_meshsets.end();mit++) {
    ierr = mField.get_msId_3dENTS_sides(*mit,true,0); CHKERRQ(ierr);
    ierr = mField.get_msId_3dENTS_split_sides(0,bit_level0,*mit,true,true,0); CHKERRQ(ierr);
  }
  EntityHandle meshset_level0;
  rval = moab.create_meshset(MESHSET_SET,meshset_level0); CHKERR_PETSC(rval);
  ierr = mField.refine_get_ents(bit_level0,meshset_level0); CHKERRQ(ierr);

  BitRefLevel bit_level1;
  bit_level1.set(1);
  ierr = mField.seed_ref_level_3D(meshset_level0,bit_level1); CHKERRQ(ierr);

  EntityHandle meshset_level1;
  rval = moab.create_meshset(MESHSET_SET,meshset_level1); CHKERR_PETSC(rval);
  ierr = mField.refine_get_ents(bit_level1,meshset_level1); CHKERRQ(ierr);

  // random mesh refinment
  EntityHandle meshset_ref_edges;
  rval = moab.create_meshset(MESHSET_SET,meshset_ref_edges); CHKERR_PETSC(rval);
  Range edges_to_refine;
  rval = moab.get_entities_by_type(meshset_level1,MBEDGE,edges_to_refine);  CHKERR_PETSC(rval);
  for(Range::iterator eit = edges_to_refine.begin();eit!=edges_to_refine.end();eit++) {
    int numb = rand() % 3;
    if(numb == 0) {
      ierr = moab.add_entities(meshset_ref_edges,&*eit,1); CHKERRQ(ierr);
    }
  }
  BitRefLevel bit_level2;
  bit_level2.set(2);
  ierr = mField.add_verices_in_the_middel_of_edges(meshset_ref_edges,bit_level2); CHKERRQ(ierr);
  ierr = mField.refine_TET(meshset_level1,bit_level2); CHKERRQ(ierr);
  ierr = mField.refine_PRISM(meshset_level1,bit_level2,0); CHKERRQ(ierr);

  ierr = PetscGetTime(&v2);CHKERRQ(ierr);
  ierr = PetscGetCPUTime(&t2);CHKERRQ(ierr);
  PetscSynchronizedPrintf(PETSC_COMM_WORLD,"SNESSolve:: Rank %d Time = %f CPU Time = %f\n",pcomm->rank(),v2-v1,t2-t1);
  PetscSynchronizedFlush(PETSC_COMM_WORLD);

  EntityHandle out_meshset_ref;
  rval = moab.create_meshset(MESHSET_SET,out_meshset_ref); CHKERR_PETSC(rval);
  ierr = mField.refine_get_finite_elements(bit_level2,out_meshset_ref); CHKERRQ(ierr);
  rval = moab.write_file(mesh_out_file_name,"VTK","",&out_meshset_ref,1); CHKERR_PETSC(rval);

  PetscFinalize();

}


