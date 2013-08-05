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

#include "moabField.hpp"
#include "moabField_Core.hpp"
#include "moabFEMethod_UpLevelStudent.hpp"

using namespace MoFEM;

ErrorCode rval;
PetscErrorCode ierr;

static char help[] = "...\n\n";

int main(int argc, char *argv[]) {

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
  char mesh_out_file_name[255];
  ierr = PetscOptionsGetString(PETSC_NULL,"-my_out_file",mesh_out_file_name,255,&flg); CHKERRQ(ierr);
  if(flg != PETSC_TRUE) {
    SETERRQ(PETSC_COMM_SELF,1,"*** ERROR -my_out_file (MESH FILE NEEDED)");
  }

  const char *option;
  option = "";//"PARALLEL=BCAST;";//;DEBUG_IO";
  rval = moab.load_file(mesh_file_name, 0, option); CHKERR_PETSC(rval); 
  ParallelComm* pcomm = ParallelComm::get_pcomm(&moab,MYPCOMM_INDEX);
  if(pcomm == NULL) pcomm =  new ParallelComm(&moab,PETSC_COMM_WORLD);

  PetscLogDouble t1,t2;
  PetscLogDouble v1,v2;
  ierr = PetscGetTime(&v1); CHKERRQ(ierr);
  ierr = PetscGetCPUTime(&t1); CHKERRQ(ierr);

  moabField_Core core(moab);
  moabField& mField = core;

  Range CubitSideSets_meshsets;
  ierr = mField.get_CubitBCType_meshsets(SideSet,CubitSideSets_meshsets); CHKERRQ(ierr);

  //ref meshset
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

  //ref level 1
  BitRefLevel bit_level1;
  bit_level1.set(1);
  ierr = mField.add_verices_in_the_middel_of_edges(meshset_level0,bit_level1); CHKERRQ(ierr);
  ierr = mField.refine_TET(meshset_level0,bit_level1); CHKERRQ(ierr);
  ierr = mField.refine_PRISM(meshset_level0,bit_level1); CHKERRQ(ierr);

  EntityHandle meshset_level1;
  rval = moab.create_meshset(MESHSET_SET,meshset_level1); CHKERR_PETSC(rval);
  ierr = mField.refine_get_ents(bit_level1,meshset_level1); CHKERRQ(ierr);

  //add fields
  ierr = mField.add_field("SPATIAL_POSITION",H1,3); CHKERRQ(ierr);
  ierr = mField.add_field("DAMAGE",L2,1); CHKERRQ(ierr);
  ierr = mField.add_field("LAMBDA",NoField,1); CHKERRQ(ierr);
  //ierr = mField.add_field("DAMAGE_INTERFACE",L2_2D,1); CHKERRQ(ierr);

  //add finite elements
  ierr = mField.add_finite_element("ELASTIC"); CHKERRQ(ierr);
  ierr = mField.modify_finite_element_add_field_row("ELASTIC","SPATIAL_POSITION"); CHKERRQ(ierr);
  ierr = mField.modify_finite_element_add_field_col("ELASTIC","SPATIAL_POSITION"); CHKERRQ(ierr);
  ierr = mField.modify_finite_element_add_field_data("ELASTIC","SPATIAL_POSITION"); CHKERRQ(ierr);
  ierr = mField.modify_finite_element_add_field_data("ELASTIC","DAMAGE"); CHKERRQ(ierr);
  ierr = mField.add_finite_element("FE_DAMAGE"); CHKERRQ(ierr);
  ierr = mField.modify_finite_element_add_field_row("FE_DAMAGE","DAMAGE"); CHKERRQ(ierr);
  ierr = mField.modify_finite_element_add_field_col("FE_DAMAGE","DAMAGE"); CHKERRQ(ierr);
  ierr = mField.modify_finite_element_add_field_data("FE_DAMAGE","DAMAGE"); CHKERRQ(ierr);
  ierr = mField.modify_finite_element_add_field_data("FE_DAMAGE","SPATIAL_POSITION"); CHKERRQ(ierr);
  ierr = mField.add_finite_element("COUPLING_WITH_DAMAGE1"); CHKERRQ(ierr);
  ierr = mField.modify_finite_element_add_field_row("COUPLING_WITH_DAMAGE1","DAMAGE"); CHKERRQ(ierr);
  ierr = mField.modify_finite_element_add_field_col("COUPLING_WITH_DAMAGE1","SPATIAL_POSITION"); CHKERRQ(ierr);
  ierr = mField.modify_finite_element_add_field_data("COUPLING_WITH_DAMAGE1","SPATIAL_POSITION"); CHKERRQ(ierr);
  ierr = mField.modify_finite_element_add_field_data("COUPLING_WITH_DAMAGE1","DAMAGE"); CHKERRQ(ierr);
  ierr = mField.add_finite_element("COUPLING_WITH_DAMAGE2"); CHKERRQ(ierr);
  ierr = mField.modify_finite_element_add_field_row("COUPLING_WITH_DAMAGE2","SPATIAL_POSITION"); CHKERRQ(ierr);
  ierr = mField.modify_finite_element_add_field_col("COUPLING_WITH_DAMAGE2","DAMAGE"); CHKERRQ(ierr);
  ierr = mField.modify_finite_element_add_field_data("COUPLING_WITH_DAMAGE2","SPATIAL_POSITION"); CHKERRQ(ierr);
  ierr = mField.modify_finite_element_add_field_data("COUPLING_WITH_DAMAGE2","DAMAGE"); CHKERRQ(ierr);
  /*ierr = mField.add_finite_element("ARC_LENGHT"); CHKERRQ(ierr);
  ierr = mField.modify_finite_element_add_field_row("ARC_LENGHT","LAMBDA"); CHKERRQ(ierr);
  ierr = mField.modify_finite_element_add_field_col("ARC_LENGHT","LAMBDA"); CHKERRQ(ierr);
  ierr = mField.modify_finite_element_add_field_col("ARC_LENGHT","SPATIAL_POSITION"); CHKERRQ(ierr);
  ierr = mField.modify_finite_element_add_field_data("ARC_LENGHT","LAMBDA"); CHKERRQ(ierr);
  ierr = mField.modify_finite_element_add_field_data("ARC_LENGHT","SPATIAL_POSITION"); CHKERRQ(ierr);
  ierr = mField.add_finite_element("LINK"); CHKERRQ(ierr);
  ierr = mField.modify_finite_element_add_field_row("LINK","SPATIAL_POSITION"); CHKERRQ(ierr);
  ierr = mField.modify_finite_element_add_field_col("LINK","SPATIAL_POSITION"); CHKERRQ(ierr);
  ierr = mField.modify_finite_element_add_field_data("LINK","SPATIAL_POSITION"); CHKERRQ(ierr);*/
  ierr = mField.add_finite_element("INTERFACE"); CHKERRQ(ierr);
  ierr = mField.modify_finite_element_add_field_row("INTERFACE","SPATIAL_POSITION"); CHKERRQ(ierr);
  //ierr = mField.modify_finite_element_add_field_row("INTERFACE","DAMAGE_INTERFACE"); CHKERRQ(ierr);
  ierr = mField.modify_finite_element_add_field_col("INTERFACE","SPATIAL_POSITION"); CHKERRQ(ierr);
  //ierr = mField.modify_finite_element_add_field_col("INTERFACE","DAMAGE_INTERFACE"); CHKERRQ(ierr);
  ierr = mField.modify_finite_element_add_field_data("INTERFACE","SPATIAL_POSITION"); CHKERRQ(ierr);
  //ierr = mField.modify_finite_element_add_field_data("INTERFACE","DAMAGE_INTERFACE"); CHKERRQ(ierr);

  //add problems 
  ierr = mField.add_problem("ELASTIC_MECHANICS"); CHKERRQ(ierr);
  ierr = mField.add_problem("DAMAGE_MECHANICS"); CHKERRQ(ierr);
  ierr = mField.add_problem("COUPLED_DAMAGE_MECHANICS"); CHKERRQ(ierr);
  //ierr = mField.add_problem("LINK_PROBLEM"); CHKERRQ(ierr);
  //ierr = mField.add_problem("ARC_LENGHT_PROBLEM"); CHKERRQ(ierr);
  ierr = mField.add_problem("ELASTIC_MECHANICS_LEVEL0"); CHKERRQ(ierr);


  //define problems and finite elements
  ierr = mField.modify_problem_add_finite_element("ELASTIC_MECHANICS","ELASTIC"); CHKERRQ(ierr);
  //ierr = mField.modify_problem_add_finite_element("ELASTIC_MECHANICS","LINK"); CHKERRQ(ierr);
  ierr = mField.modify_problem_add_finite_element("DAMAGE_MECHANICS","FE_DAMAGE"); CHKERRQ(ierr);
  ierr = mField.modify_problem_add_finite_element("COUPLED_DAMAGE_MECHANICS","ELASTIC"); CHKERRQ(ierr);
  ierr = mField.modify_problem_add_finite_element("COUPLED_DAMAGE_MECHANICS","FE_DAMAGE"); CHKERRQ(ierr);
  ierr = mField.modify_problem_add_finite_element("COUPLED_DAMAGE_MECHANICS","COUPLING_WITH_DAMAGE1"); CHKERRQ(ierr);
  ierr = mField.modify_problem_add_finite_element("COUPLED_DAMAGE_MECHANICS","COUPLING_WITH_DAMAGE2"); CHKERRQ(ierr);
  //ierr = mField.modify_problem_add_finite_element("COUPLED_DAMAGE_MECHANICS","LINK"); CHKERRQ(ierr);
  //ierr = mField.modify_problem_add_finite_element("COUPLED_DAMAGE_MECHANICS","ARC_LENGHT"); CHKERRQ(ierr);
  ierr = mField.modify_problem_add_finite_element("ELASTIC_MECHANICS_LEVEL0","ELASTIC"); CHKERRQ(ierr);
  ierr = mField.modify_problem_add_finite_element("ELASTIC_MECHANICS_LEVEL0","INTERFACE"); CHKERRQ(ierr);

  /*//special elments
  ierr = mField.modify_problem_add_finite_element("LINK_PROBLEM","LINK"); CHKERRQ(ierr);
  ierr = mField.modify_problem_add_finite_element("ARC_LENGHT_PROBLEM","ARC_LENGHT"); CHKERRQ(ierr);*/
  
  ierr = mField.modify_problem_ref_level_add_bit("ELASTIC_MECHANICS",bit_level1); CHKERRQ(ierr);
  ierr = mField.modify_problem_ref_level_add_bit("COUPLED_DAMAGE_MECHANICS",bit_level1); CHKERRQ(ierr);
  ierr = mField.modify_problem_ref_level_add_bit("ELASTIC_MECHANICS_LEVEL0",bit_level0); CHKERRQ(ierr);

  //add entities to the field meshset
  ierr = mField.add_ents_to_field_by_TETs(0,"DAMAGE"); CHKERRQ(ierr);
  ierr = mField.add_ents_to_field_by_TETs(0,"SPATIAL_POSITION"); CHKERRQ(ierr);
  //ierr = mField.add_ents_to_field_by_PRISMs(0,"DAMAGE_INTERFACE"); CHKERRQ(ierr);

  //add finite elements entities
  ierr = mField.add_ents_to_finite_element_EntType_by_bit_ref(bit_level0,"ELASTIC",MBTET); CHKERRQ(ierr);
  ierr = mField.add_ents_to_finite_element_EntType_by_bit_ref(bit_level0,"FE_DAMAGE",MBTET); CHKERRQ(ierr);
  ierr = mField.add_ents_to_finite_element_EntType_by_bit_ref(bit_level0,"COUPLING_WITH_DAMAGE1",MBTET); CHKERRQ(ierr);
  ierr = mField.add_ents_to_finite_element_EntType_by_bit_ref(bit_level0,"COUPLING_WITH_DAMAGE2",MBTET); CHKERRQ(ierr);
  ierr = mField.add_ents_to_finite_element_EntType_by_bit_ref(bit_level1,"ELASTIC",MBTET); CHKERRQ(ierr);
  ierr = mField.add_ents_to_finite_element_EntType_by_bit_ref(bit_level1,"FE_DAMAGE",MBTET); CHKERRQ(ierr);
  ierr = mField.add_ents_to_finite_element_EntType_by_bit_ref(bit_level1,"COUPLING_WITH_DAMAGE1",MBTET); CHKERRQ(ierr);
  ierr = mField.add_ents_to_finite_element_EntType_by_bit_ref(bit_level1,"COUPLING_WITH_DAMAGE2",MBTET); CHKERRQ(ierr);
  ierr = mField.add_ents_to_finite_element_EntType_by_bit_ref(bit_level0,"INTERFACE",MBPRISM); CHKERRQ(ierr);


  /*ierr = mField.add_ents_to_finite_element_by_MESHSET(Block5,"LINK"); CHKERRQ(ierr);
  EntityHandle meshset_FE_ARC_LENGHT;
  rval = moab.create_meshset(MESHSET_SET,meshset_FE_ARC_LENGHT); CHKERR_PETSC(rval);
  EntityHandle meshset_field_LAMBDA = mField.get_field_meshset("LAMBDA");
  rval = moab.add_entities(meshset_FE_ARC_LENGHT,&meshset_field_LAMBDA,1); CHKERR_PETSC(rval);
  rval = moab.add_entities(meshset_FE_ARC_LENGHT,&Block5,1); CHKERR_PETSC(rval);
  ierr = mField.seed_ref_level_MESHSET(meshset_FE_ARC_LENGHT,BitRefLevel().set()); CHKERRQ(ierr);
  ierr = mField.add_ents_to_finite_element_by_MESHSET(meshset_FE_ARC_LENGHT,"ARC_LENGHT"); CHKERRQ(ierr);*/

  //set app. order
  ierr = mField.set_field_order(0,MBTET,"SPATIAL_POSITION",4); CHKERRQ(ierr);
  ierr = mField.set_field_order(0,MBTET,"DAMAGE",0); CHKERRQ(ierr);
  ierr = mField.set_field_order(0,MBTRI,"SPATIAL_POSITION",4); CHKERRQ(ierr);
  ierr = mField.set_field_order(0,MBEDGE,"SPATIAL_POSITION",4); CHKERRQ(ierr);
  ierr = mField.set_field_order(0,MBVERTEX,"SPATIAL_POSITION",1); CHKERRQ(ierr);
  //ierr = mField.set_field_order(0,MBPRISM,"DAMAGE_INTERFACE",1); CHKERRQ(ierr);

  moabField_Core core2(moab);
  moabField& mField2 = core2;

  //build fields
  ierr = mField2.build_fields(); CHKERRQ(ierr);
  ierr = mField2.build_fields(); CHKERRQ(ierr);

  //build finite elements
  ierr = mField2.build_finite_elements(); CHKERRQ(ierr);
  ierr = mField2.build_finite_elements(); CHKERRQ(ierr);

  //build adjacencies
  ierr = mField2.build_adjacencies(bit_level0|bit_level1); CHKERRQ(ierr);
  ierr = mField2.build_adjacencies(bit_level0|bit_level1); CHKERRQ(ierr);

  //build problem
  ierr = mField2.build_problems(); CHKERRQ(ierr);

  //partition
  ierr = mField2.partition_problems("COUPLED_DAMAGE_MECHANICS"); CHKERRQ(ierr);
  ierr = mField2.partition_problems("ELASTIC_MECHANICS"); CHKERRQ(ierr);
  ierr = mField2.partition_problems("ELASTIC_MECHANICS_LEVEL0"); CHKERRQ(ierr);

  ierr = mField2.partition_finite_elements("COUPLED_DAMAGE_MECHANICS"); CHKERRQ(ierr);
  ierr = mField2.partition_finite_elements("ELASTIC_MECHANICS"); CHKERRQ(ierr);
  ierr = mField2.partition_finite_elements("ELASTIC_MECHANICS_LEVEL0"); CHKERRQ(ierr);

  ierr = mField2.partition_ghost_dofs("COUPLED_DAMAGE_MECHANICS"); CHKERRQ(ierr);
  ierr = mField2.partition_ghost_dofs("ELASTIC_MECHANICS"); CHKERRQ(ierr);
  ierr = mField2.partition_ghost_dofs("ELASTIC_MECHANICS_LEVEL0"); CHKERRQ(ierr);

  Vec Dofs_row_ELASTIC_MECHANICS;
  ierr = mField2.VecCreateGhost("ELASTIC_MECHANICS",Row,&Dofs_row_ELASTIC_MECHANICS); CHKERRQ(ierr);
  Vec Dofs_col_ELASTIC_MECHANICS;
  ierr = mField2.VecCreateGhost("ELASTIC_MECHANICS",Col,&Dofs_col_ELASTIC_MECHANICS); CHKERRQ(ierr);
  Vec Dofs_row_ELASTIC_MECHANICS_LEVEL0;
  ierr = mField2.VecCreateGhost("ELASTIC_MECHANICS_LEVEL0",Row,&Dofs_row_ELASTIC_MECHANICS_LEVEL0); CHKERRQ(ierr);
  Mat Aij_ELASTIC_MECHANICS_LEVEL0;
  ierr = mField2.MatCreateMPIAIJWithArrays("ELASTIC_MECHANICS_LEVEL0",&Aij_ELASTIC_MECHANICS_LEVEL0); CHKERRQ(ierr);

  struct MyFEMethod: public FEMethod_UpLevelStudent {

    Mat &Aij;
    Vec& rows_vec;
    MyFEMethod(Interface& _moab,Mat &_Aij,Vec& _rows_vec): FEMethod_UpLevelStudent(_moab,1),Aij(_Aij),rows_vec(_rows_vec) { }; 

    vector<double> g_NTET;

    PetscErrorCode preProcess() {
      PetscFunctionBegin;
      g_NTET.resize(4*4);
      ShapeMBTET(&g_NTET[0],G_TET_X4,G_TET_Y4,G_TET_Z4,4);
      PetscFunctionReturn(0);
    }


    PetscErrorCode operator()() {
      PetscFunctionBegin;
      ierr = OpStudentStart_TET(g_NTET); CHKERRQ(ierr);

      ostringstream ss;

      vector<DofIdx> RowGlobDofs;
      ierr = GetRowGlobalIndices("SPATIAL_POSITION",RowGlobDofs); CHKERRQ(ierr);
      vector<DofIdx> ColGlobDofs;
      ierr = GetColGlobalIndices("SPATIAL_POSITION",ColGlobDofs); CHKERRQ(ierr);

      vector<vector<DofIdx> > RowGlobEdges;
      RowGlobEdges.resize(6);
      for(int ee = 0;ee<6;ee++) {
	ierr = GetRowGlobalIndices("SPATIAL_POSITION",MBEDGE,RowGlobEdges[ee],ee); CHKERRQ(ierr);
      }
      
      vector<ublas::matrix<double> > nodeRowNMatrix;
      ierr = GetGaussRowNMatrix("SPATIAL_POSITION",nodeRowNMatrix); CHKERRQ(ierr);
      vector<ublas::matrix<double> > nodeRowDiffNMatrix;
      ierr = GetGaussRowDiffNMatrix("SPATIAL_POSITION",nodeRowDiffNMatrix); CHKERRQ(ierr);
      vector<ublas::matrix<double> > edgeRowNMatrix[6];
      vector<ublas::matrix<double> > edgeRowDiffNMatrix[6];
      for(int ee = 0;ee<6;ee++) {
	ierr = GetGaussRowNMatrix("SPATIAL_POSITION",MBEDGE,edgeRowNMatrix[ee],ee); CHKERRQ(ierr);
	ierr = GetGaussRowDiffNMatrix("SPATIAL_POSITION",MBEDGE,edgeRowDiffNMatrix[ee],ee); CHKERRQ(ierr);
      }
      vector<ublas::matrix<double> > faceRowNMatrix[4];
      vector<ublas::matrix<double> > faceRowDiffNMatrix[4];
      for(int ff = 0;ff<4;ff++) {
	ierr = GetGaussRowNMatrix("SPATIAL_POSITION",MBTRI,faceRowNMatrix[ff],ff); CHKERRQ(ierr);
	ierr = GetGaussRowDiffNMatrix("SPATIAL_POSITION",MBTRI,faceRowDiffNMatrix[ff],ff); CHKERRQ(ierr);
      }

      vector<ublas::matrix<double> > nodeBMatrix;
      ierr = MakeBMatrix3D("SPATIAL_POSITION",nodeRowDiffNMatrix,nodeBMatrix); CHKERRQ(ierr);
      vector<ublas::matrix<double> > edgeBMatrix0;
      ierr = MakeBMatrix3D("SPATIAL_POSITION",edgeRowDiffNMatrix[0],edgeBMatrix0); CHKERRQ(ierr);
  
      ublas::matrix<FieldData> NTN,edgeNTN[6],faceNTN[4];
      int g_dim = g_NTET.size()/4;
      for(int gg = 0;gg<g_dim;gg++) {
	NTN = prod(trans(nodeRowDiffNMatrix[gg]), nodeRowDiffNMatrix[gg]);
	for(int ee = 0;ee<6;ee++) {
	  edgeNTN[ee] = prod(trans(edgeRowNMatrix[ee][gg]),edgeRowNMatrix[ee][gg]);
	}
	for(int ff = 0;ff<4;ff++) {
	  faceNTN[ff] = prod(trans(faceRowNMatrix[ff][gg]),faceRowNMatrix[ff][gg]);
	}
      }


      //copy(RowGlobDofs.begin(),RowGlobDofs.end(),ostream_iterator<DofIdx>(ss,", ")); ss << endl;
      //copy(ColGlobDofs.begin(),ColGlobDofs.end(),ostream_iterator<DofIdx>(ss,", ")); ss << endl;
      //for(int ee = 0;ee<6;ee++) {
	//copy(RowGlobEdges[ee].begin(),RowGlobEdges[ee].end(),ostream_iterator<DofIdx>(ss,", ")); ss << endl;
      //}
      PetscSynchronizedPrintf(PETSC_COMM_WORLD,ss.str().c_str());

      ierr = OpStudentEnd(); CHKERRQ(ierr);
      PetscFunctionReturn(0);
    }
    PetscErrorCode postProcess() {
      PetscFunctionBegin;
      PetscFunctionReturn(0);
    }
  };

  MyFEMethod method(moab,Aij_ELASTIC_MECHANICS_LEVEL0,Dofs_row_ELASTIC_MECHANICS_LEVEL0);
  ierr = mField2.loop_finite_elements("ELASTIC_MECHANICS_LEVEL0","ELASTIC",method);  CHKERRQ(ierr);
  PetscSynchronizedFlush(PETSC_COMM_WORLD);

  ierr = VecDestroy(&Dofs_row_ELASTIC_MECHANICS); CHKERRQ(ierr);
  ierr = VecDestroy(&Dofs_col_ELASTIC_MECHANICS); CHKERRQ(ierr);
  ierr = VecDestroy(&Dofs_row_ELASTIC_MECHANICS_LEVEL0); CHKERRQ(ierr);
  ierr = MatDestroy(&Aij_ELASTIC_MECHANICS_LEVEL0); CHKERRQ(ierr);

  ierr = PetscGetTime(&v2);CHKERRQ(ierr);
  ierr = PetscGetCPUTime(&t2);CHKERRQ(ierr);

  PetscSynchronizedPrintf(PETSC_COMM_WORLD,"SNESSolve:: Rank %d Time = %f CPU Time = %f\n",pcomm->rank(),v2-v1,t2-t1);
  PetscSynchronizedFlush(PETSC_COMM_WORLD);

  EntityHandle out_meshset;
  rval = moab.create_meshset(MESHSET_SET,out_meshset); CHKERR_PETSC(rval);
  //ierr = mField2.refine_get_finite_elements(bit_level1,out_meshset); CHKERRQ(ierr);
  ierr = mField2.problem_get_FE("ELASTIC_MECHANICS_LEVEL0","ELASTIC",out_meshset); CHKERRQ(ierr);

  moab.write_file(mesh_out_file_name,"VTK","",&out_meshset,1);

  PetscFinalize();

}

