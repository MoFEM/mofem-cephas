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
  ParallelComm* pcomm = ParallelComm::get_pcomm(&moab,MYPCOMM_INDEX);
  if(pcomm == NULL) pcomm =  new ParallelComm(&moab,PETSC_COMM_WORLD);

  FieldCore core(moab);
  FieldInterface& mField = core;

  ierr = mField.seed_ref_level_3D(0,BitRefLevel().set(0)); CHKERRQ(ierr);
  vector<BitRefLevel> bit_levels;
  bit_levels.push_back(BitRefLevel().set(0));

  int ll = 1;
  for(_IT_CUBITMESHSETS_BY_BCDATA_TYPE_FOR_LOOP_(mField,SideSet|InterfaceSet,cit)) {
    ierr = PetscPrintf(PETSC_COMM_WORLD,"Insert Interface %d\n",cit->get_msId()); CHKERRQ(ierr);
    EntityHandle cubit_meshset = cit->get_meshset();
    {
      //get tet enties form back bit_level
      EntityHandle ref_level_meshset = 0;
      rval = moab.create_meshset(MESHSET_SET,ref_level_meshset); CHKERR_PETSC(rval);
      ierr = mField.refine_get_ents(bit_levels.back(),BitRefLevel().set(),MBTET,ref_level_meshset); CHKERRQ(ierr);
      ierr = mField.refine_get_ents(bit_levels.back(),BitRefLevel().set(),MBPRISM,ref_level_meshset); CHKERRQ(ierr);
      Range ref_level_tets;
      rval = moab.get_entities_by_handle(ref_level_meshset,ref_level_tets,true); CHKERR_PETSC(rval);
      //get faces and test to split
      ierr = mField.get_msId_3dENTS_sides(cubit_meshset,true,0); CHKERRQ(ierr);
      //set new bit level
      bit_levels.push_back(BitRefLevel().set(ll++));
      //split faces and edges
      ierr = mField.get_msId_3dENTS_split_sides(ref_level_meshset,bit_levels.back(),cubit_meshset,true,true,0); CHKERRQ(ierr);
      //clean meshsets
      rval = moab.delete_entities(&ref_level_meshset,1); CHKERR_PETSC(rval);
    }
    //update cubit meshsets
    for(_IT_CUBITMESHSETS_FOR_LOOP_(mField,ciit)) {
      EntityHandle cubit_meshset = ciit->meshset; 
      ierr = mField.refine_get_childern(cubit_meshset,bit_levels.back(),cubit_meshset,MBVERTEX,true); CHKERRQ(ierr);
      ierr = mField.refine_get_childern(cubit_meshset,bit_levels.back(),cubit_meshset,MBEDGE,true); CHKERRQ(ierr);
      ierr = mField.refine_get_childern(cubit_meshset,bit_levels.back(),cubit_meshset,MBTRI,true); CHKERRQ(ierr);
      ierr = mField.refine_get_childern(cubit_meshset,bit_levels.back(),cubit_meshset,MBTET,true); CHKERRQ(ierr);
    }
  }

  //add filds
  ierr = mField.add_field("H1FIELD_SCALAR",H1,1); CHKERRQ(ierr);

  //add finite elements
  ierr = mField.add_finite_element("ELEM_SCALAR"); CHKERRQ(ierr);
  ierr = mField.modify_finite_element_add_field_row("ELEM_SCALAR","H1FIELD_SCALAR"); CHKERRQ(ierr);
  ierr = mField.modify_finite_element_add_field_col("ELEM_SCALAR","H1FIELD_SCALAR"); CHKERRQ(ierr);
  ierr = mField.modify_finite_element_add_field_data("ELEM_SCALAR","H1FIELD_SCALAR"); CHKERRQ(ierr);
  //FE Interface
  ierr = mField.add_finite_element("INTERFACE"); CHKERRQ(ierr);
  ierr = mField.modify_finite_element_add_field_row("INTERFACE","H1FIELD_SCALAR"); CHKERRQ(ierr);
  ierr = mField.modify_finite_element_add_field_col("INTERFACE","H1FIELD_SCALAR"); CHKERRQ(ierr);
  ierr = mField.modify_finite_element_add_field_data("INTERFACE","H1FIELD_SCALAR"); CHKERRQ(ierr);

  //add ents to field and set app. order
  ierr = mField.add_ents_to_field_by_TETs(0,"H1FIELD_SCALAR"); CHKERRQ(ierr);
  ierr = mField.set_field_order(0,MBVERTEX,"H1FIELD_SCALAR",1); CHKERRQ(ierr);

  //add finite elements entities
  //all TETS and PRIMS are added to finite elements, for testin pruposes.
  //in some practiacl applications to save memory, you would like to add elements
  //from particular refinment level (see: mField.add_ents_to_finite_element_EntType_by_bit_ref(...)
  ierr = mField.add_ents_to_finite_element_by_TETs(0,"ELEM_SCALAR",MBTET); CHKERRQ(ierr);
  ierr = mField.add_ents_to_finite_element_by_TETs(0,"INTERFACE",MBPRISM); CHKERRQ(ierr);

  //add problems 
  //set problem for all levels, only for testing pruposes
  for(int lll = 0;lll<ll;lll++) {
    stringstream problem_name;
    problem_name << "PROBLEM_SCALAR_" << lll;
    ierr = mField.add_problem(problem_name.str()); CHKERRQ(ierr);
    //define problems and finite elements
    ierr = mField.modify_problem_add_finite_element(problem_name.str(),"ELEM_SCALAR"); CHKERRQ(ierr);
    ierr = mField.modify_problem_add_finite_element(problem_name.str(),"INTERFACE"); CHKERRQ(ierr);
  }

  //set problem level
  for(int lll = 0;lll<ll;lll++) {
    stringstream problem_name;
    problem_name << "PROBLEM_SCALAR_" << lll;
    stringstream message;
    message << "set problem problem < " << problem_name.str() << " > bit level " << bit_levels[lll] << endl;
    ierr = PetscPrintf(PETSC_COMM_WORLD,message.str().c_str()); CHKERRQ(ierr);
    ierr = mField.modify_problem_ref_level_add_bit(problem_name.str(),bit_levels[lll]); CHKERRQ(ierr);
  }

  //build fields
  ierr = mField.build_fields(); CHKERRQ(ierr);
  //build finite elements
  ierr = mField.build_finite_elements(); CHKERRQ(ierr);
  //build adjacencies
  //Its build adjacencies for all ements in databse,
  //for pratical applications consider to build adjacencies
  //only for refinemnt levels which you use for calulations
  ierr = mField.build_adjacencies(BitRefLevel().set()); CHKERRQ(ierr);
  //build problem
  ierr = mField.build_problems(); CHKERRQ(ierr);

  //partition
  for(int lll = 0;lll<ll;lll++) {
    stringstream problem_name;
    problem_name << "PROBLEM_SCALAR_" << lll;
    ierr = mField.partition_problem(problem_name.str()); CHKERRQ(ierr);
    ierr = mField.partition_finite_elements(problem_name.str()); CHKERRQ(ierr);
    ierr = mField.partition_ghost_dofs(problem_name.str()); CHKERRQ(ierr);
  }

  ofstream myfile;
  myfile.open("mesh_insert_interface.txt");

  EntityHandle out_meshset_tet;
  rval = moab.create_meshset(MESHSET_SET,out_meshset_tet); CHKERR_PETSC(rval);
  ierr = mField.refine_get_ents(bit_levels.back(),BitRefLevel().set(),MBTET,out_meshset_tet); CHKERRQ(ierr);
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
  EntityHandle out_meshset_prism;
  rval = moab.create_meshset(MESHSET_SET,out_meshset_prism); CHKERR_PETSC(rval);
  ierr = mField.refine_get_ents(bit_levels.back(),BitRefLevel().set(),MBPRISM,out_meshset_prism); CHKERRQ(ierr);
  Range prisms;
  rval = moab.get_entities_by_handle(out_meshset_prism,prisms); CHKERR_PETSC(rval);
  for(Range::iterator pit = prisms.begin();pit!=prisms.end();pit++) {
    int num_nodes; 
    const EntityHandle* conn;
    rval = moab.get_connectivity(*pit,conn,num_nodes,true); CHKERR_PETSC(rval);
    
    for(int nn = 0;nn<num_nodes;nn++) {
      cout << conn[nn] << " ";
      myfile << conn[nn] << " ";
    }
    cout << endl;
    myfile << endl;

  }
  myfile.close();

  //rval = moab.write_file(mesh_out_file_name,"VTK","",&out_meshset_tet,1); CHKERR_PETSC(rval);
  //rval = moab.write_file("out_tet.vtk","VTK","",&out_meshset_tet,1); CHKERR_PETSC(rval);
  //rval = moab.write_file("out_prism.vtk","VTK","",&out_meshset_prism,1); CHKERR_PETSC(rval);

  PetscFinalize();

}
