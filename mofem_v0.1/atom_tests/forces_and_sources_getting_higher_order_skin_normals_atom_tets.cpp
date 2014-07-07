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
#include <moab/Skinner.hpp>

#include "Projection10NodeCoordsOnField.hpp"

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
  ierr = mField.add_field("FIELD1",H1,3); CHKERRQ(ierr);

  //FE
  ierr = mField.add_finite_element("TEST_FE"); CHKERRQ(ierr);

  //Define rows/cols and element data
  ierr = mField.modify_finite_element_add_field_row("TEST_FE","FIELD1"); CHKERRQ(ierr);
  ierr = mField.modify_finite_element_add_field_col("TEST_FE","FIELD1"); CHKERRQ(ierr);
  ierr = mField.modify_finite_element_add_field_data("TEST_FE","FIELD1"); CHKERRQ(ierr);

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
  //add entities to finite element
  Range tets;
  rval = moab.get_entities_by_type(0,MBTET,tets,false); CHKERR_PETSC(rval);
  Skinner skin(&mField.get_moab());
  Range tets_skin;
  rval = skin.find_skin(tets,false,tets_skin); CHKERR(rval);
  ierr = mField.add_ents_to_finite_element_by_TRIs(tets_skin,"TEST_FE"); CHKERRQ(ierr);

  //set app. order
  //see Hierarchic Finite Element Bases on Unstructured Tetrahedral Meshes (Mark Ainsworth & Joe Coyle)
  int order = 3;
  ierr = mField.set_field_order(root_set,MBTET,"FIELD1",order); CHKERRQ(ierr);
  ierr = mField.set_field_order(root_set,MBTRI,"FIELD1",order); CHKERRQ(ierr);
  ierr = mField.set_field_order(root_set,MBEDGE,"FIELD1",order); CHKERRQ(ierr);
  ierr = mField.set_field_order(root_set,MBVERTEX,"FIELD1",1); CHKERRQ(ierr);

  /****/
  //build database
  //build field
  ierr = mField.build_fields(); CHKERRQ(ierr);
  //set FIELD1 from positions of 10 node tets
  Projection10NodeCoordsOnField ent_method(mField,"FIELD1");
  ierr = mField.loop_dofs("FIELD1",ent_method); CHKERRQ(ierr);
  //build finite elemnts
  ierr = mField.build_finiteElementsPtr(); CHKERRQ(ierr);
  //build adjacencies
  ierr = mField.build_adjacencies(bit_level0); CHKERRQ(ierr);
  //build problem
  ierr = mField.build_problems(); CHKERRQ(ierr);

  /****/
  //mesh partitioning 
  //partition
  ierr = mField.simple_partition_problem("TEST_PROBLEM"); CHKERRQ(ierr);
  ierr = mField.partition_finiteElementsPtr("TEST_PROBLEM"); CHKERRQ(ierr);
  //what are ghost nodes, see Petsc Manual
  ierr = mField.partition_ghost_dofs("TEST_PROBLEM"); CHKERRQ(ierr);

  struct ForcesAndSurcesCore_TestFE: public ForcesAndSurcesCore {

    ErrorCode rval;
    PetscErrorCode ierr;

    ublas::matrix<FieldData> nOrmals_at_GaussPt;
    ublas::matrix<FieldData> tAngent1_at_GaussPt;
    ublas::matrix<FieldData> tAngent2_at_GaussPt;
    DataForcesAndSurcesCore data;
    OpGetNormals op;

    typedef tee_device<ostream, ofstream> TeeDevice;
    typedef stream<TeeDevice> TeeStream;
    ofstream ofs;
    TeeDevice my_tee; 
    TeeStream my_split;

    ForcesAndSurcesCore_TestFE(FieldInterface &_mField): 
      ForcesAndSurcesCore(_mField),data(MBTRI),
      op(nOrmals_at_GaussPt,tAngent1_at_GaussPt,tAngent2_at_GaussPt),
      ofs("forces_and_sources_getting_higher_order_skin_normals_atom.txt"),
      my_tee(cout,ofs),my_split(my_tee) {};

    PetscErrorCode preProcess() {
      PetscFunctionBegin;
      PetscFunctionReturn(0);
    }



    PetscErrorCode operator()() {
      PetscFunctionBegin;

      ierr = getEdgesSense(data); CHKERRQ(ierr);
      ierr = getEdgesOrder(data); CHKERRQ(ierr);
      ierr = getTrisOrder(data); CHKERRQ(ierr);
      ierr = getRowNodesIndices(data,"FIELD1"); CHKERRQ(ierr);
      ierr = getEdgesRowIndices(data,"FIELD1"); CHKERRQ(ierr);
      ierr = getTrisRowIndices(data,"FIELD1"); CHKERRQ(ierr);
      ierr = getNodesFieldData(data,"FIELD1"); CHKERRQ(ierr);
      ierr = getEdgesFieldData(data,"FIELD1"); CHKERRQ(ierr);
      ierr = getTrisFieldData(data,"FIELD1"); CHKERRQ(ierr);

      ierr = shapeTRIFunctions_H1(data,G_TRI_X4,G_TRI_Y4,4); CHKERRQ(ierr);

      nOrmals_at_GaussPt.resize(4,3);
      tAngent1_at_GaussPt.resize(4,3);
      tAngent2_at_GaussPt.resize(4,3);

      try {
	ierr = op.op(data); CHKERRQ(ierr);
	ierr = op.calculateNormals(); CHKERRQ(ierr);
      } catch (exception& ex) {
	ostringstream ss;
	ss << "thorw in method: " << ex.what() << " at line " << __LINE__ << " in file " << __FILE__ << endl;
	SETERRQ(PETSC_COMM_SELF,1,ss.str().c_str());
      }

      my_split.precision(3);
      my_split << "normals: " << nOrmals_at_GaussPt << endl;
      my_split << "tangent1: " << tAngent1_at_GaussPt << endl;
      my_split << "tangent2: " << tAngent2_at_GaussPt << endl;

      PetscFunctionReturn(0);
    }

    PetscErrorCode postProcess() {
      PetscFunctionBegin;
      PetscFunctionReturn(0);
    }


  };

  ForcesAndSurcesCore_TestFE fe1(mField);
  ierr = mField.loop_finiteElementsPtr("TEST_PROBLEM","TEST_FE",fe1);  CHKERRQ(ierr);

  ierr = PetscFinalize(); CHKERRQ(ierr);

  return 0;

}


