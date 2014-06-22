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
#include "FEMethod_UpLevelStudent.hpp"
#include "PotentialFlowFEMethod.hpp"
#include "SurfacePressure.hpp"

#include <petscksp.h>
#include "Projection10NodeCoordsOnField.hpp"
#include "PostProcDisplacementAndStrainOnRefindedMesh.hpp"

using namespace MoFEM;

ErrorCode rval;
PetscErrorCode ierr;

static char help[] = "...\n\n";


int main(int argc, char *argv[]) {

  try {

  PetscInitialize(&argc,&argv,(char *)0,help);

  Core mb_instance;
  Interface& moab = mb_instance;
  int rank;
  MPI_Comm_rank(PETSC_COMM_WORLD,&rank);
  /*if(rank==0) {
    EntityHandle dummy_meshset;
    rval = moab.create_meshset(MESHSET_SET,dummy_meshset); CHKERR_PETSC(rval);
  }*/

  PetscBool flg = PETSC_TRUE;
  char mesh_file_name[255];
  ierr = PetscOptionsGetString(PETSC_NULL,"-my_file",mesh_file_name,255,&flg); CHKERRQ(ierr);
  if(flg != PETSC_TRUE) {
    SETERRQ(PETSC_COMM_SELF,1,"*** ERROR -my_file (MESH FILE NEEDED)");
  }
  PetscInt order;
  ierr = PetscOptionsGetInt(PETSC_NULL,"-my_order",&order,&flg); CHKERRQ(ierr);
  if(flg != PETSC_TRUE) {
    order = 1;
  }


  const char *option;
  option = "";//"PARALLEL=BCAST;";//;DEBUG_IO";
  rval = moab.load_file(mesh_file_name, 0, option); CHKERR(rval); 
  ParallelComm* pcomm = ParallelComm::get_pcomm(&moab,MYPCOMM_INDEX);
  if(pcomm == NULL) pcomm =  new ParallelComm(&moab,PETSC_COMM_WORLD);

  FieldCore core(moab);
  FieldInterface& mField = core;

  //add filds
  ierr = mField.add_field("POTENTIAL_FIELD",H1,1); CHKERRQ(ierr);
  ierr = mField.add_field("MESH_NODE_POSITIONS",H1,3); CHKERRQ(ierr);

  //add finite elements
  ierr = mField.add_finite_element("POTENTIAL_ELEM"); CHKERRQ(ierr);
  ierr = mField.modify_finite_element_add_field_row("POTENTIAL_ELEM","POTENTIAL_FIELD"); CHKERRQ(ierr);
  ierr = mField.modify_finite_element_add_field_col("POTENTIAL_ELEM","POTENTIAL_FIELD"); CHKERRQ(ierr);
  ierr = mField.modify_finite_element_add_field_data("POTENTIAL_ELEM","POTENTIAL_FIELD"); CHKERRQ(ierr);
  ierr = mField.modify_finite_element_add_field_data("POTENTIAL_ELEM","MESH_NODE_POSITIONS"); CHKERRQ(ierr);

  //add problems 
  ierr = mField.add_problem("POTENTIAL_PROBLEM"); CHKERRQ(ierr);

  //define problems and finite elements
  ierr = mField.modify_problem_add_finite_element("POTENTIAL_PROBLEM","POTENTIAL_ELEM"); CHKERRQ(ierr);

  BitRefLevel bit_level0;
  bit_level0.set(0);
  ierr = mField.seed_ref_level_3D(0,bit_level0); CHKERRQ(ierr);

  for(_IT_CUBITMESHSETS_BY_BCDATA_TYPE_FOR_LOOP_(mField,BlockSet|UnknownCubitName,it)) {
    if(it->get_Cubit_name() == "PotentialFlow") {
 
      //add ents to field and set app. order
      ierr = mField.add_ents_to_field_by_TETs(0,"POTENTIAL_FIELD"); CHKERRQ(ierr);
      ierr = mField.add_ents_to_field_by_TETs(0,"MESH_NODE_POSITIONS"); CHKERRQ(ierr);

      //add finite elements entities
      ierr = mField.add_ents_to_finite_element_by_TETs(it->meshset,"POTENTIAL_ELEM",true); CHKERRQ(ierr);

    }
  }

  ierr = mField.set_field_order(0,MBVERTEX,"POTENTIAL_FIELD",1); CHKERRQ(ierr);
  ierr = mField.set_field_order(0,MBEDGE,"POTENTIAL_FIELD",order); CHKERRQ(ierr);
  ierr = mField.set_field_order(0,MBTRI,"POTENTIAL_FIELD",order); CHKERRQ(ierr);
  ierr = mField.set_field_order(0,MBTET,"POTENTIAL_FIELD",order); CHKERRQ(ierr);
  //
  ierr = mField.set_field_order(0,MBTET,"MESH_NODE_POSITIONS",2); CHKERRQ(ierr);
  ierr = mField.set_field_order(0,MBTRI,"MESH_NODE_POSITIONS",2); CHKERRQ(ierr);
  ierr = mField.set_field_order(0,MBEDGE,"MESH_NODE_POSITIONS",2); CHKERRQ(ierr);
  ierr = mField.set_field_order(0,MBVERTEX,"MESH_NODE_POSITIONS",1); CHKERRQ(ierr);

  //flux boundary conditions
  ierr = MetaNeummanForces::addNeumannFluxBCElements(mField,"POTENTIAL_PROBLEM","POTENTIAL_FIELD"); CHKERRQ(ierr);

  //set problem level
  ierr = mField.modify_problem_ref_level_add_bit("POTENTIAL_PROBLEM",bit_level0); CHKERRQ(ierr);

  //build fields
  ierr = mField.build_fields(); CHKERRQ(ierr);
  //build finite elements
  ierr = mField.build_finite_elements(); CHKERRQ(ierr);
  //build adjacencies
  ierr = mField.build_adjacencies(bit_level0); CHKERRQ(ierr);
  //build problem
  ierr = mField.build_problems(); CHKERRQ(ierr);

  //partition problems
  ierr = mField.partition_problem("POTENTIAL_PROBLEM"); CHKERRQ(ierr);
  ierr = mField.partition_finite_elements("POTENTIAL_PROBLEM"); CHKERRQ(ierr);
  ierr = mField.partition_ghost_dofs("POTENTIAL_PROBLEM"); CHKERRQ(ierr);

  //print bcs
  ierr = mField.print_cubit_displacement_set(); CHKERRQ(ierr);
  ierr = mField.print_cubit_pressure_set(); CHKERRQ(ierr);

  //create matrices and vectors
  Vec F,D;
  ierr = mField.VecCreateGhost("POTENTIAL_PROBLEM",Row,&F); CHKERRQ(ierr);
  ierr = mField.VecCreateGhost("POTENTIAL_PROBLEM",Col,&D); CHKERRQ(ierr);
  Mat A;
  ierr = mField.MatCreateMPIAIJWithArrays("POTENTIAL_PROBLEM",&A); CHKERRQ(ierr);

  Projection10NodeCoordsOnField ent_method_material(mField,"MESH_NODE_POSITIONS");
  ierr = mField.loop_dofs("MESH_NODE_POSITIONS",ent_method_material); CHKERRQ(ierr);

  //get nodes and other entities to fix
  Range fix_nodes;
  for(_IT_CUBITMESHSETS_BY_BCDATA_TYPE_FOR_LOOP_(mField,NodeSet|UnknownCubitName,it)) {
    std::size_t zeroPressureFound=it->get_Cubit_name().find("ZeroPressure");
    if (zeroPressureFound==std::string::npos) continue;
    rval = moab.get_entities_by_type(it->meshset,MBVERTEX,fix_nodes,true); CHKERR_PETSC(rval);
    Range edges;
    rval = moab.get_entities_by_type(it->meshset,MBEDGE,edges,true); CHKERR_PETSC(rval);
    Range tris;
    rval = moab.get_entities_by_type(it->meshset,MBTRI,tris,true); CHKERR_PETSC(rval);
    Range adj;
    rval = moab.get_connectivity(tris,adj,true); CHKERR_PETSC(rval);
    fix_nodes.insert(adj.begin(),adj.end());
    rval = moab.get_connectivity(edges,adj,true); CHKERR_PETSC(rval);
    fix_nodes.insert(adj.begin(),adj.end());
    rval = moab.get_adjacencies(tris,1,false,edges,Interface::UNION); CHKERR_PETSC(rval);
  }
  FixMaterialPoints fix_dofs(mField,"POTENTIAL_FIELD",A,D,F,fix_nodes);
  //initialize data structure
  ierr = mField.problem_basic_method_preProcess("POTENTIAL_PROBLEM",fix_dofs); CHKERRQ(ierr);

  //neuman flux bc elements
  boost::ptr_map<string,NeummanForcesSurface> neumann_forces;
  ierr = MetaNeummanForces::setNeumannFluxFiniteElementOperators(mField,neumann_forces,F,"POTENTIAL_FIELD"); CHKERRQ(ierr);
  boost::ptr_map<string,NeummanForcesSurface>::iterator mit = neumann_forces.begin();
  for(;mit!=neumann_forces.end();mit++) {
    ierr = mField.loop_finite_elements("POTENTIAL_PROBLEM",mit->first,mit->second->getLoopFe()); CHKERRQ(ierr);
  }
  //evaluate laplacian in body
  LaplacianElem elem(mField,A,F);
  ierr = MatZeroEntries(A); CHKERRQ(ierr);
  ierr = mField.loop_finite_elements("POTENTIAL_PROBLEM","POTENTIAL_ELEM",elem);  CHKERRQ(ierr);
  
  //post proces fix boundary conditiond
  ierr = mField.problem_basic_method_postProcess("POTENTIAL_PROBLEM",fix_dofs); CHKERRQ(ierr);

  ierr = MatAssemblyBegin(A,MAT_FINAL_ASSEMBLY); CHKERRQ(ierr);
  ierr = MatAssemblyEnd(A,MAT_FINAL_ASSEMBLY); CHKERRQ(ierr);
  ierr = VecAssemblyBegin(F); CHKERRQ(ierr);
  ierr = VecAssemblyEnd(F); CHKERRQ(ierr);

  //Matrix View
  /*{
    MatView(A,PETSC_VIEWER_DRAW_WORLD);//PETSC_VIEWER_STDOUT_WORLD);
    std::string wait;
    std::cin >> wait;
  }*/
  //ierr = VecView(F,PETSC_VIEWER_STDOUT_WORLD); CHKERRQ(ierr);

  //Solver
  KSP solver;
  ierr = KSPCreate(PETSC_COMM_WORLD,&solver); CHKERRQ(ierr);
  ierr = KSPSetOperators(solver,A,A,SAME_NONZERO_PATTERN); CHKERRQ(ierr);
  ierr = KSPSetFromOptions(solver); CHKERRQ(ierr);
  ierr = KSPSetUp(solver); CHKERRQ(ierr);

  ierr = KSPSolve(solver,F,D); CHKERRQ(ierr);
  ierr = VecGhostUpdateBegin(D,INSERT_VALUES,SCATTER_FORWARD); CHKERRQ(ierr);
  ierr = VecGhostUpdateEnd(D,INSERT_VALUES,SCATTER_FORWARD); CHKERRQ(ierr);
  ierr = mField.set_global_VecCreateGhost("POTENTIAL_PROBLEM",Row,D,INSERT_VALUES,SCATTER_REVERSE); CHKERRQ(ierr);

  ierr = KSPDestroy(&solver); CHKERRQ(ierr);
  ierr = VecDestroy(&F); CHKERRQ(ierr);
  ierr = VecDestroy(&D); CHKERRQ(ierr);
  ierr = MatDestroy(&A); CHKERRQ(ierr);

  /*Tag th_phi;
  double def_val = 0;
  rval = moab.tag_get_handle("PHI",1,MB_TYPE_DOUBLE,th_phi,MB_TAG_CREAT|MB_TAG_SPARSE,&def_val); CHKERR_PETSC(rval);
  for(_IT_GET_DOFS_FIELD_BY_NAME_FOR_LOOP_(mField,"POTENTIAL_FIELD",dof)) {
    EntityHandle ent = dof->get_ent();
    double val = dof->get_FieldData();
    rval = moab.tag_set_data(th_phi,&ent,1,&val); CHKERR_PETSC(rval);
  }*/

  ProjectionFieldOn10NodeTet ent_method_phi_on_10nodeTet(mField,"POTENTIAL_FIELD",true,false,"PHI");
  ierr = mField.loop_dofs("POTENTIAL_FIELD",ent_method_phi_on_10nodeTet); CHKERRQ(ierr);
  ent_method_phi_on_10nodeTet.set_nodes = false;
  ierr = mField.loop_dofs("POTENTIAL_FIELD",ent_method_phi_on_10nodeTet); CHKERRQ(ierr);

  if(pcomm->rank()==0) {
    rval = moab.write_file("solution.h5m"); CHKERR_PETSC(rval);
  }

  if(pcomm->rank()==0) {
    EntityHandle out_meshset;
    rval = moab.create_meshset(MESHSET_SET,out_meshset); CHKERR_PETSC(rval);
    ierr = mField.problem_get_FE("POTENTIAL_PROBLEM","POTENTIAL_ELEM",out_meshset); CHKERRQ(ierr);
    rval = moab.write_file("out.vtk","VTK","",&out_meshset,1); CHKERR_PETSC(rval);
    rval = moab.delete_entities(&out_meshset,1); CHKERR_PETSC(rval);
  }

  PostProcScalarFieldsAndGradientOnRefMesh fe_post_proc_method(moab,"POTENTIAL_FIELD");
  ierr = mField.loop_finite_elements("POTENTIAL_PROBLEM","POTENTIAL_ELEM",fe_post_proc_method);  CHKERRQ(ierr);

  if(pcomm->rank()==0) {
    rval = fe_post_proc_method.moab_post_proc.write_file("out_post_proc.vtk","VTK",""); CHKERR_PETSC(rval);
  }

  PetscFinalize();

  } catch (const char* msg) {
    SETERRQ(PETSC_COMM_SELF,1,msg);
  } catch (const std::exception& ex) {
    ostringstream ss;
    ss << "thorw in method: " << ex.what() << endl;
    SETERRQ(PETSC_COMM_SELF,1,ss.str().c_str());
  }

}
