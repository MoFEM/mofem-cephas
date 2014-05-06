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
#include "ForcesAndSurcesCore.hpp"

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
  int rank;
  MPI_Comm_rank(PETSC_COMM_WORLD,&rank);

  PetscBool flg = PETSC_TRUE;
  char mesh_file_name[255];
  ierr = PetscOptionsGetString(PETSC_NULL,"-my_file",mesh_file_name,255,&flg); CHKERRQ(ierr);
  if(flg != PETSC_TRUE) {
    SETERRQ(PETSC_COMM_SELF,1,"*** ERROR -my_file (MESH FILE NEEDED)");
  }

  //Create MoFEM (Joseph) database
  FieldCore core(moab);
  FieldInterface& mField = core;

  ParallelComm* pcomm = ParallelComm::get_pcomm(&moab,MYPCOMM_INDEX);
  if(pcomm == NULL) pcomm =  new ParallelComm(&moab,PETSC_COMM_WORLD);

  const char *option;
  option = "";//"PARALLEL=BCAST;";//;DEBUG_IO";
  BARRIER_RANK_START(pcomm) 
  rval = moab.load_file(mesh_file_name, 0, option); CHKERR_PETSC(rval); 
  BARRIER_RANK_END(pcomm) 

  //set entitities bit level
  BitRefLevel bit_level0;
  bit_level0.set(0);
  EntityHandle meshset_level0;
  rval = moab.create_meshset(MESHSET_SET,meshset_level0); CHKERR_PETSC(rval);
  ierr = mField.seed_ref_level_3D(0,bit_level0); CHKERRQ(ierr);

  //Fields
  ierr = mField.add_field("FIELD1",H1,1); CHKERRQ(ierr);
  ierr = mField.add_field("FIELD2",H1,3); CHKERRQ(ierr);

  //FE
  ierr = mField.add_finite_element("TEST_FE"); CHKERRQ(ierr);

  //Define rows/cols and element data
  ierr = mField.modify_finite_element_add_field_row("TEST_FE","FIELD1"); CHKERRQ(ierr);
  ierr = mField.modify_finite_element_add_field_col("TEST_FE","FIELD2"); CHKERRQ(ierr);
  ierr = mField.modify_finite_element_add_field_data("TEST_FE","FIELD1"); CHKERRQ(ierr);
  ierr = mField.modify_finite_element_add_field_data("TEST_FE","FIELD2"); CHKERRQ(ierr);

  //Problem
  ierr = mField.add_problem("TEST_PROBLEM"); CHKERRQ(ierr);

  //set finite elements for problem
  ierr = mField.modify_problem_add_finite_element("TEST_PROBLEM","TEST_FE"); CHKERRQ(ierr);
  //set refinment level for problem
  ierr = mField.modify_problem_ref_level_add_bit("TEST_PROBLEM",bit_level0); CHKERRQ(ierr);


  //meshset consisting all entities in mesh
  EntityHandle root_set = moab.get_root_set(); 
  //add entities to field
  ierr = mField.add_ents_to_field_by_TETs(root_set,"FIELD1"); CHKERRQ(ierr);
  ierr = mField.add_ents_to_field_by_TETs(root_set,"FIELD2"); CHKERRQ(ierr);
  //add entities to finite element
  ierr = mField.add_ents_to_finite_element_by_TETs(root_set,"TEST_FE"); CHKERRQ(ierr);


  //set app. order
  //see Hierarchic Finite Element Bases on Unstructured Tetrahedral Meshes (Mark Ainsworth & Joe Coyle)
  int order = 5;
  ierr = mField.set_field_order(root_set,MBTET,"FIELD1",order); CHKERRQ(ierr);
  ierr = mField.set_field_order(root_set,MBTRI,"FIELD1",order); CHKERRQ(ierr);
  ierr = mField.set_field_order(root_set,MBEDGE,"FIELD1",order); CHKERRQ(ierr);
  ierr = mField.set_field_order(root_set,MBVERTEX,"FIELD1",1); CHKERRQ(ierr);
  ierr = mField.set_field_order(root_set,MBTET,"FIELD2",order); CHKERRQ(ierr);
  ierr = mField.set_field_order(root_set,MBTRI,"FIELD2",order); CHKERRQ(ierr);
  ierr = mField.set_field_order(root_set,MBEDGE,"FIELD2",order); CHKERRQ(ierr);
  ierr = mField.set_field_order(root_set,MBVERTEX,"FIELD2",1); CHKERRQ(ierr);

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

  struct ForcesAndSurcesCore_TestFE: public ForcesAndSurcesCore {

    ErrorCode rval;
    PetscErrorCode ierr;

    typedef tee_device<ostream, ofstream> TeeDevice;
    typedef stream<TeeDevice> TeeStream;

    ofstream ofs;
    TeeDevice my_tee; 
    TeeStream my_split;


    ForcesAndSurcesCore_TestFE(FieldInterface &_mField): 
      ForcesAndSurcesCore(_mField), 
      ofs("forces_and_surces_getting_orders_indices_atom_test.txt"),
      my_tee(cout, ofs),my_split(my_tee) {};

    PetscErrorCode preProcess() {
      PetscFunctionBegin;
      PetscFunctionReturn(0);
    }

    PetscErrorCode operator()() {
      PetscFunctionBegin;
      my_split << "\n\nNEXT ELEM\n\n";
      ierr = getEdgesSense(); CHKERRQ(ierr);
      my_split << "edgesSense: " << edgesSense << endl;
      ierr = getFacesSense(); CHKERRQ(ierr);
      my_split << "facesSense: " << facesSense << endl;
      ierr = getEdgesOrder("FIELD1"); CHKERRQ(ierr);
      my_split << "FIELD1 edgesOrder: " << edgesOrder << endl;
      ierr = getEdgesOrder("FIELD2"); CHKERRQ(ierr);
      my_split << "FIELD2 edgesOrder: " << edgesOrder << endl;
      ierr = getFacesOrder("FIELD1"); CHKERRQ(ierr);
      my_split << "FIELD1 facesOrder: " << facesOrder << endl;
      ierr = getFacesOrder("FIELD2"); CHKERRQ(ierr);
      my_split << "FIELD2 facesOrder: " << facesOrder << endl;
      ierr = getOrderVolume("FIELD1"); CHKERRQ(ierr);
      my_split << "FIELD1 volumeOrder: " << volumeOrder << endl;
      ierr = getOrderVolume("FIELD2"); CHKERRQ(ierr);
      my_split << "FIELD2 volumeOrder: " << volumeOrder << endl;
      ierr = getRowNodesIndices("FIELD1"); CHKERRQ(ierr);
      my_split << "FIELD1 rowNodesIndices: " << rowNodesIndices << endl;
      ierr = getColNodesIndices("FIELD2"); CHKERRQ(ierr);
      my_split << "FIELD2 colNodesIndices: " << colNodesIndices << endl;
      ierr = getEdgeRowIndices("FIELD1"); CHKERRQ(ierr);
      for(int ee = 0;ee<6;ee++) {
	my_split << "FIELD1 " << ee << " " << rowEdgesIndcies[ee] << endl;
      }
      ierr = getEdgeColIndices("FIELD2"); CHKERRQ(ierr);
      for(int ee = 0;ee<6;ee++) {
	my_split << "FIELD2 " << ee << " " << colEdgesIndcies[ee] << endl;
      }
      ierr = getFacesRowIndices("FIELD1"); CHKERRQ(ierr);
      for(int ff = 0;ff<4;ff++) {
	my_split << "FIELD2 " << ff << " " <<  rowFacesIndcies[ff] << endl;
      }
      ierr = getFacesColIndices("FIELD2"); CHKERRQ(ierr);
      for(int ff = 0;ff<4;ff++) {
	my_split << "FIELD2 " << ff << " " <<  colFacesIndcies[ff] << endl;
      }
      ierr = getTetRowIndices("FIELD1"); CHKERRQ(ierr);
      my_split << "FIELD1 " << rowTetIndcies << endl;
      ierr = getTetColIndices("FIELD2"); CHKERRQ(ierr);
      my_split << "FIELD1 " << colTetIndcies << endl;

      ierr = getFaceNodes(); CHKERRQ(ierr);
      my_split << "facesNodes " << facesNodes << endl;

      ierr = shapeTETFunctions_H1(G_TET_X4,G_TET_Y4,G_TET_Z4,4); CHKERRQ(ierr);
      my_split.precision(2);
      my_split << "Nnodes_H1 " << Nnodes_H1 << endl;
      my_split << "diffNodes_H1 " << diffNnodes_H1 << endl;
      for(int ee = 0;ee<6;ee++) {
	my_split << "Nedges_H1 ee = " << ee << " " << Nedges_H1[ee] << endl;
      }
      for(int ff = 0;ff<4;ff++) {
	my_split << "Nedges_H1 ff = " << ff << " " << Nfaces_H1[ff] << endl;
      }
      my_split << "Nvolume_H1 " << Nvolume_H1 << endl;

      ierr = shapeTETFunctions_L2(G_TET_X4,G_TET_Y4,G_TET_Z4,4); CHKERRQ(ierr);
      my_split << "Nnodes_L2 " << Nnodes_H1 << endl;
      my_split << "diffNodes_L2 " << diffNnodes_H1 << endl;
      my_split << "Nvolume_L2 " << Nvolume_H1 << endl;

      PetscFunctionReturn(0);
    }

    PetscErrorCode postProcess() {
      PetscFunctionBegin;

      my_split.close();

      PetscFunctionReturn(0);
    }

  };

  ForcesAndSurcesCore_TestFE fe1(mField);
  ierr = mField.loop_finite_elements("TEST_PROBLEM","TEST_FE",fe1);  CHKERRQ(ierr);

  ierr = PetscFinalize(); CHKERRQ(ierr);

  return 0;

}


