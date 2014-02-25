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

  Core mb_instance;
  Interface& moab = mb_instance;

  const char *option;
  option = "";//"PARALLEL=BCAST";//;DEBUG_IO";
  rval = moab.load_file(mesh_file_name, 0, option); CHKERR(rval); 

  FieldCore core(moab);
  FieldInterface& mField = core;

  BitRefLevel bit_level0;
  bit_level0.set(0);
  ierr = mField.seed_ref_level_3D(0,bit_level0); CHKERRQ(ierr);

  BitRefLevel bit_level1;
  bit_level1.set(1);

  EntityHandle meshset_level0;
  rval = moab.create_meshset(MESHSET_SET,meshset_level0); CHKERR_PETSC(rval);
  ierr = mField.get_entities_by_type_and_ref_level(bit_level0,BitRefLevel().set(),meshset_level0); CHKERRQ(ierr);

  // random mesh refinment
  EntityHandle meshset_ref_edges;
  rval = moab.create_meshset(MESHSET_SET,meshset_ref_edges); CHKERR_PETSC(rval);
  Range edges_to_refine;
  rval = moab.get_entities_by_type(meshset_level0,MBEDGE,edges_to_refine);  CHKERR_PETSC(rval);
  int ii = 0;
  for(Range::iterator eit = edges_to_refine.begin();
    eit!=edges_to_refine.end();eit++,ii++) {
    int numb = ii % 2;
    if(numb == 0) {
      ierr = moab.add_entities(meshset_ref_edges,&*eit,1); CHKERRQ(ierr);
    }
  }
  ierr = mField.add_verices_in_the_middel_of_edges(meshset_ref_edges,bit_level1); CHKERRQ(ierr);
  ierr = mField.refine_TET(meshset_level0,bit_level1); CHKERRQ(ierr);

  ierr = mField.shift_right_bit_ref(1); CHKERRQ(ierr);

  ofstream myfile;
  myfile.open("mesh_refine.txt");

  EntityHandle out_meshset_tet;
  rval = moab.create_meshset(MESHSET_SET,out_meshset_tet); CHKERR_PETSC(rval);
  ierr = mField.get_entities_by_type_and_ref_level(bit_level0,BitRefLevel().set(),MBTET,out_meshset_tet); CHKERRQ(ierr);
  Range tets;
  rval = moab.get_entities_by_handle(out_meshset_tet,tets); CHKERR_PETSC(rval);
  for(Range::iterator tit = tets.begin();tit!=tets.end();tit++) {
    int num_nodes; 
    const EntityHandle* conn;
    rval = moab.get_connectivity(*tit,conn,num_nodes,true); CHKERR_PETSC(rval);
    
    for(int nn = 0;nn<num_nodes;nn++) {
      cout << conn[nn] << " ";
      myfile << conn[nn] << " ";
    }
    cout << endl;
    myfile << endl;

  }

  myfile.close();

  //rval = moab.write_file("out.vtk","VTK","",&out_meshset_tet,1); CHKERR_PETSC(rval);

  PetscFinalize();

}
