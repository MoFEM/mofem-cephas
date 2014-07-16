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

#include "SurfacePressure.hpp"
#include <boost/ptr_container/ptr_map.hpp>

#include <boost/iostreams/tee.hpp>
#include <boost/iostreams/stream.hpp>
#include <fstream>
#include <iostream>

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
  ierr = mField.add_field("DISPLACEMENT",H1,3); CHKERRQ(ierr);
  ierr = mField.add_field("MESH_NODE_POSITIONS",H1,3); CHKERRQ(ierr);

  //Problem
  ierr = mField.add_problem("TEST_PROBLEM"); CHKERRQ(ierr);
  //set refinment level for problem
  ierr = mField.modify_problem_ref_level_add_bit("TEST_PROBLEM",bit_level0); CHKERRQ(ierr);

  //meshset consisting all entities in mesh
  EntityHandle root_set = moab.get_root_set(); 
  //add entities to field
  ierr = mField.add_ents_to_field_by_TETs(root_set,"DISPLACEMENT"); CHKERRQ(ierr);
  ierr = mField.add_ents_to_field_by_TETs(root_set,"MESH_NODE_POSITIONS"); CHKERRQ(ierr);

  for(_IT_CUBITMESHSETS_BY_BCDATA_TYPE_FOR_LOOP_(mField,NODESET|FORCESET,it)) {

    ostringstream fe_name;
    fe_name << "FORCE_FE_" << it->get_msId();
    ierr = mField.add_finite_element(fe_name.str()); CHKERRQ(ierr);
    ierr = mField.modify_finite_element_add_field_row(fe_name.str(),"DISPLACEMENT"); CHKERRQ(ierr);
    ierr = mField.modify_finite_element_add_field_col(fe_name.str(),"DISPLACEMENT"); CHKERRQ(ierr);
    ierr = mField.modify_finite_element_add_field_data(fe_name.str(),"DISPLACEMENT"); CHKERRQ(ierr);
    ierr = mField.modify_problem_add_finite_element("TEST_PROBLEM",fe_name.str()); CHKERRQ(ierr);

    Range tris;
    rval = moab.get_entities_by_type(it->meshset,MBTRI,tris,true); CHKERR_PETSC(rval);
    ierr = mField.add_ents_to_finite_element_by_TRIs(tris,fe_name.str()); CHKERRQ(ierr);

  }

  for(_IT_CUBITMESHSETS_BY_BCDATA_TYPE_FOR_LOOP_(mField,SIDESET|PRESSURESET,it)) {

    ostringstream fe_name;
    fe_name << "PRESSURE_FE_" << it->get_msId();
    ierr = mField.add_finite_element(fe_name.str()); CHKERRQ(ierr);
    ierr = mField.modify_finite_element_add_field_row(fe_name.str(),"DISPLACEMENT"); CHKERRQ(ierr);
    ierr = mField.modify_finite_element_add_field_col(fe_name.str(),"DISPLACEMENT"); CHKERRQ(ierr);
    ierr = mField.modify_finite_element_add_field_data(fe_name.str(),"DISPLACEMENT"); CHKERRQ(ierr);
    ierr = mField.modify_finite_element_add_field_data(fe_name.str(),"MESH_NODE_POSITIONS"); CHKERRQ(ierr);
    ierr = mField.modify_problem_add_finite_element("TEST_PROBLEM",fe_name.str()); CHKERRQ(ierr);

    Range tris;
    rval = moab.get_entities_by_type(it->meshset,MBTRI,tris,true); CHKERR_PETSC(rval);
    ierr = mField.add_ents_to_finite_element_by_TRIs(tris,fe_name.str()); CHKERRQ(ierr);

  }

  //set app. order
  //see Hierarchic Finite Element Bases on Unstructured Tetrahedral Meshes (Mark Ainsworth & Joe Coyle)
  int order = 2;
  ierr = mField.set_field_order(root_set,MBTET,"DISPLACEMENT",order); CHKERRQ(ierr);
  ierr = mField.set_field_order(root_set,MBTRI,"DISPLACEMENT",order); CHKERRQ(ierr);
  ierr = mField.set_field_order(root_set,MBEDGE,"DISPLACEMENT",order); CHKERRQ(ierr);
  ierr = mField.set_field_order(root_set,MBVERTEX,"DISPLACEMENT",1); CHKERRQ(ierr);
  ierr = mField.set_field_order(root_set,MBTET,"MESH_NODE_POSITIONS",2); CHKERRQ(ierr);
  ierr = mField.set_field_order(root_set,MBTRI,"MESH_NODE_POSITIONS",2); CHKERRQ(ierr);
  ierr = mField.set_field_order(root_set,MBEDGE,"MESH_NODE_POSITIONS",2); CHKERRQ(ierr);
  ierr = mField.set_field_order(root_set,MBVERTEX,"MESH_NODE_POSITIONS",1); CHKERRQ(ierr);

  /****/
  //build database
  //build field
  ierr = mField.build_fields(); CHKERRQ(ierr);
  //set FIELD1 from positions of 10 node tets
  Projection10NodeCoordsOnField ent_method(mField,"MESH_NODE_POSITIONS");
  ierr = mField.loop_dofs("MESH_NODE_POSITIONS",ent_method); CHKERRQ(ierr);

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

  Vec F;
  ierr = mField.VecCreateGhost("TEST_PROBLEM",ROW,&F); CHKERRQ(ierr);

  typedef tee_device<ostream, ofstream> TeeDevice;
  typedef stream<TeeDevice> TeeStream;
  ostringstream txt_name;
  txt_name << "forces_and_sources_" << mesh_file_name << ".txt";
  ofstream ofs(txt_name.str().c_str());
  TeeDevice my_tee(cout, ofs); 
  TeeStream my_split(my_tee);

  ierr = VecZeroEntries(F); CHKERRQ(ierr);
  boost::ptr_map<string,NeummanForcesSurface> neumann_forces;
  for(_IT_CUBITMESHSETS_BY_BCDATA_TYPE_FOR_LOOP_(mField,NODESET|FORCESET,it)) {
    ostringstream fe_name;
    fe_name << "FORCE_FE_" << it->get_msId();
    string fe_name_str = fe_name.str();
    neumann_forces.insert(fe_name_str,new NeummanForcesSurface(mField));
    neumann_forces.at(fe_name_str).addForce("DISPLACEMENT",F,it->get_msId()); CHKERRQ(ierr);
    ForceCubitBcData data;
    ierr = it->get_cubit_bc_data_structure(data); CHKERRQ(ierr);
    my_split << *it << endl;
    my_split << data << endl;
  }
  for(_IT_CUBITMESHSETS_BY_BCDATA_TYPE_FOR_LOOP_(mField,SIDESET|PRESSURESET,it)) {
    ostringstream fe_name;
    fe_name << "PRESSURE_FE_" << it->get_msId();
    string fe_name_str = fe_name.str();
    neumann_forces.insert(fe_name_str,new NeummanForcesSurface(mField));
    neumann_forces.at(fe_name_str).addPreassure("DISPLACEMENT",F,it->get_msId()); CHKERRQ(ierr);
    PressureCubitBcData data;
    ierr = it->get_cubit_bc_data_structure(data); CHKERRQ(ierr);
    my_split << *it << endl;
    my_split << data << endl;
  }
  boost::ptr_map<string,NeummanForcesSurface>::iterator mit = neumann_forces.begin();
  for(;mit!=neumann_forces.end();mit++) {
    ierr = mField.loop_finite_elements("TEST_PROBLEM",mit->first,mit->second->getLoopFe()); CHKERRQ(ierr);
  }
  ierr = VecAssemblyBegin(F); CHKERRQ(ierr);
  ierr = VecAssemblyEnd(F); CHKERRQ(ierr);

  ierr = mField.set_global_VecCreateGhost("TEST_PROBLEM",ROW,F,INSERT_VALUES,SCATTER_REVERSE); CHKERRQ(ierr);

  const double eps = 1e-4;

  const MoFEMProblem *problemPtr;
  ierr = mField.get_problem("TEST_PROBLEM",&problemPtr); CHKERRQ(ierr);
  for(_IT_NUMEREDDOFMOFEMENTITY_ROW_FOR_LOOP_(problemPtr,dit)) {

    my_split.precision(3);
    my_split.setf(std::ios::fixed);
    double val = fabs(dit->get_FieldData())<eps ? 0.0 : dit->get_FieldData();
    my_split << dit->get_petsc_gloabl_dof_idx() << " " << val << endl;

  }

  double sum = 0;
  ierr = VecSum(F,&sum); CHKERRQ(ierr);
  sum = fabs(sum)<eps ? 0.0 : sum;
  my_split << endl << "Sum : " << setprecision(3) << sum << endl;

  ierr = VecDestroy(&F); CHKERRQ(ierr);

  ierr = PetscFinalize(); CHKERRQ(ierr);

  return 0;

}


