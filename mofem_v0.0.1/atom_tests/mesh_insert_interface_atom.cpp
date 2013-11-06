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

static char help[] = "teting interface inserting algorithm\n\n";

int main(int argc, char *argv[]) {

  PetscInitialize(&argc,&argv,(char *)0,help);

  PetscBool flg = PETSC_TRUE;
  char mesh_file_name[255];
  ierr = PetscOptionsGetString(PETSC_NULL,"-my_file",mesh_file_name,255,&flg); CHKERRQ(ierr);
  if(flg != PETSC_TRUE) {
    SETERRQ(PETSC_COMM_SELF,1,"*** ERROR -my_file (MESH FILE NEEDED)");
  }
  /*char mesh_out_file_name[255];
  ierr = PetscOptionsGetString(PETSC_NULL,"-my_out_file",mesh_out_file_name,255,&flg); CHKERRQ(ierr);
  if(flg != PETSC_TRUE) {
    SETERRQ(PETSC_COMM_SELF,1,"*** ERROR -my_out_file (MESH FILE NEEDED)");
  }**/

  Core mb_instance;
  Interface& moab = mb_instance;
  const char *option;
  option = "";//"PARALLEL=BCAST";//;DEBUG_IO";
  rval = moab.load_file(mesh_file_name, 0, option); CHKERR(rval); 
  //ParallelComm* pcomm = ParallelComm::get_pcomm(&moab,0);

  FieldCore core(moab);
  FieldInterface& mField = core;

  ierr = mField.seed_ref_level_3D(0,0); CHKERRQ(ierr);
  vector<BitRefLevel> bit_levels;
  int ll = 0;
  for(_IT_CUBITMESHSETS_BY_BCDATA_TYPE_FOR_LOOP_(mField,SideSet|InterfaceSet,cit)) {
    ierr = PetscPrintf(PETSC_COMM_WORLD,"Insert Interface %d\n",cit->get_msId()); CHKERRQ(ierr);
    EntityHandle meshset = cit->get_meshset();
    bit_levels.push_back(BitRefLevel().set(ll++));
    ierr = mField.get_msId_3dENTS_sides(meshset,true,0); CHKERRQ(ierr);
    ierr = mField.get_msId_3dENTS_split_sides(0,bit_levels.back(),meshset,true,true,0); CHKERRQ(ierr);
    for(_IT_CUBITMESHSETS_FOR_LOOP_(mField,ciit)) {
      EntityHandle cubit_meshset = ciit->meshset; 
      ierr = mField.refine_get_childern(cubit_meshset,bit_levels.back(),cubit_meshset,MBVERTEX,true); CHKERRQ(ierr);
      ierr = mField.refine_get_childern(cubit_meshset,bit_levels.back(),cubit_meshset,MBEDGE,true); CHKERRQ(ierr);
      ierr = mField.refine_get_childern(cubit_meshset,bit_levels.back(),cubit_meshset,MBTRI,true); CHKERRQ(ierr);
      ierr = mField.refine_get_childern(cubit_meshset,bit_levels.back(),cubit_meshset,MBTET,true); CHKERRQ(ierr);
    }
  }


  ofstream myfile;
  myfile.open("mesh_insert_interface.txt");

  EntityHandle out_meshset;
  rval = moab.create_meshset(MESHSET_SET,out_meshset); CHKERR_PETSC(rval);
  ierr = mField.refine_get_ents(bit_levels.back(),BitRefLevel().set(),MBTET,out_meshset); CHKERRQ(ierr);
  Range tets;
  rval = moab.get_entities_by_handle(out_meshset,tets); CHKERR_PETSC(rval);
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

  //rval = moab.write_file(mesh_out_file_name,"VTK","",&out_meshset,1); CHKERR_PETSC(rval);

  PetscFinalize();

}
