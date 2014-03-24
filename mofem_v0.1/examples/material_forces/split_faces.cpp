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


#include "ConfigurationalFractureMechanics.hpp"
#include "FieldCore.hpp"
#include "FaceSplittinfTool.hpp"

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

  PetscLogDouble t1,t2;
  PetscLogDouble v1,v2;
  ierr = PetscTime(&v1); CHKERRQ(ierr);
  ierr = PetscGetCPUTime(&t1); CHKERRQ(ierr);

  FieldCore core(moab);
  FieldInterface& mField = core;

  ConfigurationalFractureMechanics conf_prob(mField);

  ierr = conf_prob.set_material_fire_wall(mField); CHKERRQ(ierr);

  //ref meshset ref level 0
  Tag th_my_ref_level;
  rval = mField.get_moab().tag_get_handle("_MY_REFINMENT_LEVEL",th_my_ref_level); CHKERR_PETSC(rval);
  const EntityHandle root_meshset = mField.get_moab().get_root_set();
  BitRefLevel *ptr_bit_level0;
  rval = mField.get_moab().tag_get_by_ptr(th_my_ref_level,&root_meshset,1,(const void**)&ptr_bit_level0); CHKERR_PETSC(rval);
  BitRefLevel& bit_level0 = *ptr_bit_level0;

  ierr = mField.build_fields(); CHKERRQ(ierr);
 
  FaceSplittingTools face_splitting(mField);

  BitRefLevel bit_last_ref = BitRefLevel().set(face_splitting.meshRefineBitLevels.back());
  Range tets_on_last_ref_level;
  ierr = mField.get_entities_by_type_and_ref_level(bit_last_ref,BitRefLevel().set(),MBTET,tets_on_last_ref_level); CHKERRQ(ierr);
  ierr = mField.seed_finite_elements(tets_on_last_ref_level); CHKERRQ(ierr);

  ierr = face_splitting.buildKDTreeForCrackSurface(bit_level0); CHKERRQ(ierr);

  ierr = face_splitting.meshRefine(); CHKERRQ(ierr);
  bit_last_ref = BitRefLevel().set(face_splitting.meshRefineBitLevels.back());

  ierr = face_splitting.initBitLevelData(bit_last_ref);  CHKERRQ(ierr);
  ierr = face_splitting.calculateDistanceFromCrackSurface();  CHKERRQ(ierr);
  ierr = face_splitting.getOpositeForntEdges(true); CHKERRQ(ierr);
  ierr = face_splitting.getCrackSurfaceCorssingEdges(true); CHKERRQ(ierr);

  ierr = face_splitting.getCrackFrontTets(true); CHKERRQ(ierr);
  ierr = face_splitting.chopTetsUntilNonOneLeftOnlyCrackSurfaceFaces(true,2); CHKERRQ(ierr);

  if(pcomm->rank()==0) {

    EntityHandle out_meshset;
    rval = moab.create_meshset(MESHSET_SET,out_meshset); CHKERR_PETSC(rval);
    ierr = mField.get_entities_by_type_and_ref_level(bit_last_ref,BitRefLevel().set(),MBTET,out_meshset); CHKERRQ(ierr);
    rval = moab.write_file("out.vtk","VTK","",&out_meshset,1); CHKERR_PETSC(rval);
    rval = moab.delete_entities(&out_meshset,1); CHKERR_PETSC(rval);
    //
    rval = moab.write_file("opositeFrontEdges.vtk","VTK","",&face_splitting.opositeFrontEdges,1); CHKERR_PETSC(rval);
    rval = moab.write_file("crackSurfaceCrossingEdges.vtk","VTK","",&face_splitting.crackSurfaceCrossingEdges,1); CHKERR_PETSC(rval);
    rval = moab.write_file("crackFrontTests.vtk","VTK","",&face_splitting.crackFrontTests,1); CHKERR_PETSC(rval);
    rval = moab.write_file("chopTetsFaces.vtk","VTK","",&face_splitting.chopTetsFaces,1); CHKERR_PETSC(rval);

  }

  ierr = face_splitting.addNewSurfaceFaces_to_Cubit_msId200(); CHKERRQ(ierr);
  ierr = face_splitting.splitFaces(); CHKERRQ(ierr);
  BitRefLevel bit_cat_level = BitRefLevel().set(face_splitting.meshIntefaceBitLevels.back());
  bit_level0 = bit_cat_level;

  PetscInt order;
  flg = PETSC_TRUE;
  ierr = PetscOptionsGetInt(PETSC_NULL,"-my_order",&order,&flg); CHKERRQ(ierr);
  if(flg != PETSC_TRUE) {
    order = 1;
  }
  
  ierr = face_splitting.addcrackFront_to_Cubit201(); CHKERRQ(ierr);

  ierr = mField.update_field_meshset_by_entities_children(bit_cat_level,0); CHKERRQ(ierr);
  ierr = mField.set_field_order_by_entity_type_and_bit_ref(bit_cat_level,BitRefLevel().set(),MBVERTEX,"MESH_NODE_DISPLACEMENT",1); CHKERRQ(ierr);
  ierr = mField.set_field_order_by_entity_type_and_bit_ref(bit_cat_level,BitRefLevel().set(),MBVERTEX,"MESH_NODE_POSITIONS",1); CHKERRQ(ierr);
  ierr = mField.set_field_order_by_entity_type_and_bit_ref(bit_cat_level,BitRefLevel().set(),MBVERTEX,"SPATIAL_POSITION",1); CHKERRQ(ierr);
  ierr = mField.set_field_order_by_entity_type_and_bit_ref(bit_cat_level,BitRefLevel().set(),MBVERTEX,"SPATIAL_DISPLACEMENT",1); CHKERRQ(ierr);
  ierr = mField.set_field_order_by_entity_type_and_bit_ref(bit_cat_level,BitRefLevel().set(),MBEDGE,"SPATIAL_POSITION",order); CHKERRQ(ierr);
  ierr = mField.set_field_order_by_entity_type_and_bit_ref(bit_cat_level,BitRefLevel().set(),MBEDGE,"SPATIAL_DISPLACEMENT",order); CHKERRQ(ierr);
  ierr = mField.set_field_order_by_entity_type_and_bit_ref(bit_cat_level,BitRefLevel().set(),MBTRI,"SPATIAL_POSITION",order); CHKERRQ(ierr);
  ierr = mField.set_field_order_by_entity_type_and_bit_ref(bit_cat_level,BitRefLevel().set(),MBTRI,"SPATIAL_DISPLACEMENT",order); CHKERRQ(ierr);
  ierr = mField.set_field_order_by_entity_type_and_bit_ref(bit_cat_level,BitRefLevel().set(),MBTET,"SPATIAL_POSITION",order); CHKERRQ(ierr);
  ierr = mField.set_field_order_by_entity_type_and_bit_ref(bit_cat_level,BitRefLevel().set(),MBTET,"SPATIAL_DISPLACEMENT",order); CHKERRQ(ierr);

  ierr = mField.set_field_order_by_entity_type_and_bit_ref(bit_cat_level,BitRefLevel().set(),MBVERTEX,"LAMBDA_SURFACE",1); CHKERRQ(ierr);
  ierr = mField.set_field_order_by_entity_type_and_bit_ref(bit_cat_level,BitRefLevel().set(),MBVERTEX,"LAMBDA_CRACK_SURFACE",1); CHKERRQ(ierr);

  ierr = mField.set_field_order_by_entity_type_and_bit_ref(bit_cat_level,BitRefLevel().set(),MBVERTEX,"LAMBDA_CRACK_TANGENT_CONSTRAIN",1); CHKERRQ(ierr);
  ierr = mField.set_field_order_by_entity_type_and_bit_ref(bit_cat_level,BitRefLevel().set(),MBVERTEX,"LAMBDA_CRACKFRONT_AREA",1); CHKERRQ(ierr);

  ierr = mField.update_finite_element_meshset_by_entities_children(bit_cat_level); CHKERRQ(ierr);

  ierr = mField.build_fields(); CHKERRQ(ierr);
  ierr = mField.build_finite_elements(); CHKERRQ(ierr);

  if(pcomm->rank()==0) {

    EntityHandle out_meshset;
    rval = moab.create_meshset(MESHSET_SET,out_meshset); CHKERR_PETSC(rval);
    ierr = mField.get_entities_by_type_and_ref_level(bit_cat_level,BitRefLevel().set(),MBTET,out_meshset); CHKERRQ(ierr);
    rval = moab.write_file("cat_out.vtk","VTK","",&out_meshset,1); CHKERR_PETSC(rval);
    rval = moab.delete_entities(&out_meshset,1); CHKERR_PETSC(rval);

  }

  if(pcomm->rank()==0) {

    EntityHandle meshset200;
    ierr = mField.get_Cubit_msId_meshset(200,SideSet,meshset200); CHKERRQ(ierr);
    rval = moab.write_file("cat_CrackFrontSurface.vtk","VTK","",&meshset200,1); CHKERR_PETSC(rval);
    EntityHandle meshset201;
    ierr = mField.get_Cubit_msId_meshset(201,SideSet,meshset201); CHKERRQ(ierr);
    rval = moab.write_file("cat_CrackFrontEdges.vtk","VTK","",&meshset201,1); CHKERR_PETSC(rval);

  }

  if(pcomm->rank()==0) {
    ierr = PetscPrintf(PETSC_COMM_WORLD,"Save results in out_split_faces.h5m\n"); CHKERRQ(ierr);
    rval = moab.write_file("out_split_faces.h5m"); CHKERR_PETSC(rval);
  }


  ierr = PetscTime(&v2);CHKERRQ(ierr);
  ierr = PetscGetCPUTime(&t2);CHKERRQ(ierr);
  PetscSynchronizedPrintf(PETSC_COMM_WORLD,"Total Rank %d Time = %f CPU Time = %f\n",pcomm->rank(),v2-v1,t2-t1);

  PetscFinalize();

  } catch (const char* msg) {
    SETERRQ(PETSC_COMM_SELF,1,msg);
  }

  return 0;
}
