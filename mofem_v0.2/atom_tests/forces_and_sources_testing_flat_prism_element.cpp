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


#include <MoFEM.hpp>
#include <Projection10NodeCoordsOnField.hpp>

#include <moab/Skinner.hpp>

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

  moab::Core mb_instance;
  Interface& moab = mb_instance;
  ParallelComm* pcomm = ParallelComm::get_pcomm(&moab,MYPCOMM_INDEX);
  if(pcomm == NULL) pcomm =  new ParallelComm(&moab,PETSC_COMM_WORLD);

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
  MoFEM::Core core(moab);
  FieldInterface& m_field = core;
  PrismInterface& interface = core;

  //set entitities bit level
  ierr = m_field.seed_ref_level_3D(0,BitRefLevel().set(0)); CHKERRQ(ierr);
  vector<BitRefLevel> bit_levels;
  bit_levels.push_back(BitRefLevel().set(0));

  int ll = 1;
  //for(_IT_CUBITMESHSETS_BY_BCDATA_TYPE_FOR_LOOP_(m_field,SIDESET|INTERFACESET,cit)) {
  for(_IT_CUBITMESHSETS_BY_SET_TYPE_FOR_LOOP_(m_field,SIDESET,cit)) {
    ierr = PetscPrintf(PETSC_COMM_WORLD,"Insert Interface %d\n",cit->get_msId()); CHKERRQ(ierr);
    EntityHandle cubit_meshset = cit->get_meshset();
    {
      //get tet enties form back bit_level
      EntityHandle ref_level_meshset = 0;
      rval = moab.create_meshset(MESHSET_SET,ref_level_meshset); CHKERR_PETSC(rval);
      ierr = m_field.get_entities_by_type_and_ref_level(bit_levels.back(),BitRefLevel().set(),MBTET,ref_level_meshset); CHKERRQ(ierr);
      ierr = m_field.get_entities_by_type_and_ref_level(bit_levels.back(),BitRefLevel().set(),MBPRISM,ref_level_meshset); CHKERRQ(ierr);
      Range ref_level_tets;
      rval = moab.get_entities_by_handle(ref_level_meshset,ref_level_tets,true); CHKERR_PETSC(rval);
      //get faces and test to split
      ierr = interface.get_msId_3dENTS_sides(cubit_meshset,bit_levels.back(),true,0); CHKERRQ(ierr);
      //set new bit level
      bit_levels.push_back(BitRefLevel().set(ll++));
      //split faces and 
      ierr = interface.get_msId_3dENTS_split_sides(ref_level_meshset,bit_levels.back(),cubit_meshset,true,true,0); CHKERRQ(ierr);
      //clean meshsets
      rval = moab.delete_entities(&ref_level_meshset,1); CHKERR_PETSC(rval);
    }
    //update cubit meshsets
    for(_IT_CUBITMESHSETS_FOR_LOOP_(m_field,ciit)) {
      EntityHandle cubit_meshset = ciit->meshset; 
      ierr = m_field.update_meshset_by_entities_children(cubit_meshset,bit_levels.back(),cubit_meshset,MBVERTEX,true); CHKERRQ(ierr);
      ierr = m_field.update_meshset_by_entities_children(cubit_meshset,bit_levels.back(),cubit_meshset,MBEDGE,true); CHKERRQ(ierr);
      ierr = m_field.update_meshset_by_entities_children(cubit_meshset,bit_levels.back(),cubit_meshset,MBTRI,true); CHKERRQ(ierr);
      ierr = m_field.update_meshset_by_entities_children(cubit_meshset,bit_levels.back(),cubit_meshset,MBTET,true); CHKERRQ(ierr);
    }
  }


  //Fields
  ierr = m_field.add_field("FIELD1",H1,3); CHKERRQ(ierr);
  ierr = m_field.add_field("MESH_NODE_POSITIONS",H1,3); CHKERRQ(ierr);

  //FE
  ierr = m_field.add_finite_element("TEST_FE"); CHKERRQ(ierr);

  //Define rows/cols and element data
  ierr = m_field.modify_finite_element_add_field_row("TEST_FE","FIELD1"); CHKERRQ(ierr);
  ierr = m_field.modify_finite_element_add_field_col("TEST_FE","FIELD1"); CHKERRQ(ierr);
  ierr = m_field.modify_finite_element_add_field_data("TEST_FE","FIELD1"); CHKERRQ(ierr);
  ierr = m_field.modify_finite_element_add_field_data("TEST_FE","MESH_NODE_POSITIONS"); CHKERRQ(ierr);

  //Problem
  ierr = m_field.add_problem("TEST_PROBLEM"); CHKERRQ(ierr);

  //set finite elements for problem
  ierr = m_field.modify_problem_add_finite_element("TEST_PROBLEM","TEST_FE"); CHKERRQ(ierr);
  //set refinment level for problem
  ierr = m_field.modify_problem_ref_level_add_bit("TEST_PROBLEM",bit_levels.back()); CHKERRQ(ierr);

  //meshset consisting all entities in mesh
  EntityHandle root_set = moab.get_root_set(); 
  //add entities to field
  ierr = m_field.add_ents_to_field_by_TETs(root_set,"FIELD1"); CHKERRQ(ierr);
  ierr = m_field.add_ents_to_field_by_TETs(root_set,"MESH_NODE_POSITIONS"); CHKERRQ(ierr);
  //add entities to finite element
  ierr = m_field.add_ents_to_finite_element_by_PRISMs(root_set,"TEST_FE",10); CHKERRQ(ierr);

  //set app. order
  //see Hierarchic Finite Element Bases on Unstructured Tetrahedral Meshes (Mark Ainsworth & Joe Coyle)
  int order = 3;
  ierr = m_field.set_field_order(root_set,MBTET,"FIELD1",order); CHKERRQ(ierr);
  ierr = m_field.set_field_order(root_set,MBTRI,"FIELD1",order); CHKERRQ(ierr);
  ierr = m_field.set_field_order(root_set,MBEDGE,"FIELD1",order); CHKERRQ(ierr);
  ierr = m_field.set_field_order(root_set,MBVERTEX,"FIELD1",1); CHKERRQ(ierr);

  ierr = m_field.set_field_order(root_set,MBTET,"MESH_NODE_POSITIONS",2); CHKERRQ(ierr);
  ierr = m_field.set_field_order(root_set,MBTRI,"MESH_NODE_POSITIONS",2); CHKERRQ(ierr);
  ierr = m_field.set_field_order(root_set,MBEDGE,"MESH_NODE_POSITIONS",2); CHKERRQ(ierr);
  ierr = m_field.set_field_order(root_set,MBVERTEX,"MESH_NODE_POSITIONS",1); CHKERRQ(ierr);

  /****/
  //build database
  //build field
  ierr = m_field.build_fields(); CHKERRQ(ierr);
  //set FIELD1 from positions of 10 node tets
  Projection10NodeCoordsOnField ent_method_field1(m_field,"FIELD1");
  ierr = m_field.loop_dofs("FIELD1",ent_method_field1); CHKERRQ(ierr);
  Projection10NodeCoordsOnField ent_method_mesh_positions(m_field,"MESH_NODE_POSITIONS");
  ierr = m_field.loop_dofs("MESH_NODE_POSITIONS",ent_method_mesh_positions); CHKERRQ(ierr);

  //build finite elemnts
  ierr = m_field.build_finite_elements(); CHKERRQ(ierr);
  //build adjacencies
  ierr = m_field.build_adjacencies(bit_levels.back()); CHKERRQ(ierr);
  //build problem
  ierr = m_field.build_problems(); CHKERRQ(ierr);

  /****/
  //mesh partitioning 
  //partition
  ierr = m_field.simple_partition_problem("TEST_PROBLEM"); CHKERRQ(ierr);
  ierr = m_field.partition_finite_elements("TEST_PROBLEM"); CHKERRQ(ierr);
  //what are ghost nodes, see Petsc Manual
  ierr = m_field.partition_ghost_dofs("TEST_PROBLEM"); CHKERRQ(ierr);

  FlatPrismElementForcesAndSurcesCore fe1(m_field);

  typedef tee_device<ostream, ofstream> TeeDevice;
  typedef stream<TeeDevice> TeeStream;

  ofstream ofs("forces_and_sources_testing_flat_prism_element.txt");
  TeeDevice my_tee(cout, ofs); 
  TeeStream my_split(my_tee);

  struct MyOp: public FlatPrismElementForcesAndSurcesCore::UserDataOperator {

    TeeStream &my_split;
    MyOp(TeeStream &_my_split):
      FlatPrismElementForcesAndSurcesCore::UserDataOperator("FIELD1","FIELD1"),
      my_split(_my_split) {}

    PetscErrorCode doWork(
      int side,
      EntityType type,
      DataForcesAndSurcesCore::EntData &data) {
      PetscFunctionBegin;

      if(data.getFieldData().empty()) PetscFunctionReturn(0);

      const double eps = 1e-4;
      for(
	ublas::unbounded_array<double>::iterator it = getNormal().data().begin();
	it!=getNormal().data().end();it++) {
	*it = fabs(*it)<eps ? 0.0 : *it;
      }
      for(
	ublas::unbounded_array<double>::iterator it = getNormals_at_GaussPtF3().data().begin();
	it!=getNormals_at_GaussPtF3().data().end();it++) {
	*it = fabs(*it)<eps ? 0.0 : *it;
      }
      for(
	ublas::unbounded_array<double>::iterator it = getTangent1_at_GaussPtF3().data().begin();
	it!=getTangent1_at_GaussPtF3().data().end();it++) {
	*it = fabs(*it)<eps ? 0.0 : *it;
      }
      for(
	ublas::unbounded_array<double>::iterator it = getTangent2_at_GaussPtF3().data().begin();
	it!=getTangent2_at_GaussPtF3().data().end();it++) {
	*it = fabs(*it)<eps ? 0.0 : *it;
      }
      for(
	ublas::unbounded_array<double>::iterator it = getNormals_at_GaussPtF4().data().begin();
	it!=getNormals_at_GaussPtF4().data().end();it++) {
	*it = fabs(*it)<eps ? 0.0 : *it;
      }
      for(
	ublas::unbounded_array<double>::iterator it = getTangent1_at_GaussPtF4().data().begin();
	it!=getTangent1_at_GaussPtF4().data().end();it++) {
	*it = fabs(*it)<eps ? 0.0 : *it;
      }
      for(
	ublas::unbounded_array<double>::iterator it = getTangent2_at_GaussPtF4().data().begin();
	it!=getTangent2_at_GaussPtF4().data().end();it++) {
	*it = fabs(*it)<eps ? 0.0 : *it;
      }

      my_split << "NH1" << endl;
      my_split << "side: " << side << " type: " << type << endl;
      my_split << data << endl;
      my_split << setprecision(3) << getCoords() << endl;
      my_split << setprecision(3) << getCoordsAtGaussPts() << endl;
      my_split << setprecision(3) << getArea() << endl;
      my_split << setprecision(3) << "nornal " << getNormal() << endl;
      my_split << setprecision(3) << "normal at Gauss pt F3 " << getNormals_at_GaussPtF3() << endl;
      my_split << setprecision(3) << getTangent1_at_GaussPtF3() << endl;
      my_split << setprecision(3) << getTangent2_at_GaussPtF3() << endl;
      my_split << setprecision(3) << "normal at Gauss pt F4 " << getNormals_at_GaussPtF4() << endl;
      my_split << setprecision(3) << getTangent1_at_GaussPtF4() << endl;
      my_split << setprecision(3) << getTangent2_at_GaussPtF4() << endl;
      PetscFunctionReturn(0);
    }

    PetscErrorCode doWork(
      int row_side,int col_side,
      EntityType row_type,EntityType col_type,
      DataForcesAndSurcesCore::EntData &row_data,
      DataForcesAndSurcesCore::EntData &col_data) {
      PetscFunctionBegin;

      if(row_data.getFieldData().empty()) PetscFunctionReturn(0);

      my_split << "NH1NH1" << endl;
      my_split << "row side: " << row_side << " row_type: " << row_type << endl;
      my_split << row_data << endl;
      my_split << "NH1NH1" << endl;
      my_split << "col side: " << col_side << " col_type: " << col_type << endl;
      my_split << row_data << endl;

      PetscFunctionReturn(0);
    }

  };

  fe1.get_op_to_do_Rhs().push_back(new MyOp(my_split));
  fe1.get_op_to_do_Lhs().push_back(new MyOp(my_split));

  ierr = m_field.loop_finite_elements("TEST_PROBLEM","TEST_FE",fe1);  CHKERRQ(ierr);


  ierr = PetscFinalize(); CHKERRQ(ierr);

  return 0;

}


