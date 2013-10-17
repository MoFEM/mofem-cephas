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
#include "cholesky.hpp"
#include <petscksp.h>

//PetscErrorCode  PetscGetTime(PetscLogDouble *t);

#include <petscsys.h> 
#include <petsctime.h>

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
  rval = moab.load_file(mesh_file_name, 0, option); CHKERR(rval); 
  ParallelComm* pcomm = ParallelComm::get_pcomm(&moab,MYPCOMM_INDEX);
  if(pcomm == NULL) pcomm =  new ParallelComm(&moab,PETSC_COMM_WORLD);

  PetscLogDouble t1,t2;
  PetscLogDouble v1,v2;
  ierr = PetscGetTime(&v1); CHKERRQ(ierr);
  ierr = PetscGetCPUTime(&t1); CHKERRQ(ierr);

  FieldCore core(moab);
  FieldInterface& mField = core;

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
  ierr = mField.refine_get_ents(bit_level0,BitRefLevel().set(),meshset_level0); CHKERRQ(ierr);

  //ref level 1
  BitRefLevel bit_level1;
  bit_level1.set(1);
  ierr = mField.add_verices_in_the_middel_of_edges(meshset_level0,bit_level1); CHKERRQ(ierr);
  ierr = mField.refine_TET(meshset_level0,bit_level1); CHKERRQ(ierr);
  ierr = mField.refine_PRISM(meshset_level0,bit_level1); CHKERRQ(ierr);

  EntityHandle meshset_level1;
  rval = moab.create_meshset(MESHSET_SET,meshset_level1); CHKERR_PETSC(rval);
  ierr = mField.refine_get_ents(bit_level1,BitRefLevel().set(),meshset_level1); CHKERRQ(ierr);

  //add fields
  ierr = mField.add_field("H1FIELD",H1,1); CHKERRQ(ierr);
  ierr = mField.add_field("H1FIELD_L2",L2,1); CHKERRQ(ierr);

  //add finite elements
  ierr = mField.add_finite_element("FEAPPROX"); CHKERRQ(ierr);
  ierr = mField.modify_finite_element_add_field_row("FEAPPROX","H1FIELD"); CHKERRQ(ierr);
  ierr = mField.modify_finite_element_add_field_col("FEAPPROX","H1FIELD"); CHKERRQ(ierr);
  ierr = mField.modify_finite_element_add_field_data("FEAPPROX","H1FIELD"); CHKERRQ(ierr);
  ierr = mField.add_finite_element("FEAPPROX_L2"); CHKERRQ(ierr);
  ierr = mField.modify_finite_element_add_field_row("FEAPPROX_L2","H1FIELD_L2"); CHKERRQ(ierr);
  ierr = mField.modify_finite_element_add_field_col("FEAPPROX_L2","H1FIELD_L2"); CHKERRQ(ierr);
  ierr = mField.modify_finite_element_add_field_data("FEAPPROX_L2","H1FIELD"); CHKERRQ(ierr);
  ierr = mField.modify_finite_element_add_field_data("FEAPPROX_L2","H1FIELD_L2"); CHKERRQ(ierr);

  //add problems 
  ierr = mField.add_problem("PROBLEM_APPROXIMATION"); CHKERRQ(ierr);
  ierr = mField.add_problem("PROBLEM_APPROXIMATION_REF"); CHKERRQ(ierr);

  //define problems and finite elements
  ierr = mField.modify_problem_add_finite_element("PROBLEM_APPROXIMATION","FEAPPROX"); CHKERRQ(ierr);
  ierr = mField.modify_problem_add_finite_element("PROBLEM_APPROXIMATION_REF","FEAPPROX_L2"); CHKERRQ(ierr);

  //set level
  ierr = mField.modify_problem_ref_level_add_bit("PROBLEM_APPROXIMATION",bit_level0); CHKERRQ(ierr);
  ierr = mField.modify_problem_ref_level_add_bit("PROBLEM_APPROXIMATION_REF",bit_level1); CHKERRQ(ierr);

  //add entities to the field meshset
  ierr = mField.add_ents_to_field_by_TETs(0,"H1FIELD"); CHKERRQ(ierr);
  ierr = mField.add_ents_to_field_by_TETs(0,"H1FIELD_L2"); CHKERRQ(ierr);

  //add finite elements entities
  ierr = mField.add_ents_to_finite_element_EntType_by_bit_ref(bit_level0,"FEAPPROX",MBTET); CHKERRQ(ierr);
  ierr = mField.add_ents_to_finite_element_EntType_by_bit_ref(bit_level1,"FEAPPROX",MBTET); CHKERRQ(ierr);
  ierr = mField.add_ents_to_finite_element_EntType_by_bit_ref(bit_level1,"FEAPPROX_L2",MBTET); CHKERRQ(ierr);

  //set app. order
  ierr = mField.set_field_order(0,MBTET,"H1FIELD",3); CHKERRQ(ierr);
  ierr = mField.set_field_order(0,MBTRI,"H1FIELD",3); CHKERRQ(ierr);
  ierr = mField.set_field_order(0,MBEDGE,"H1FIELD",3); CHKERRQ(ierr);
  ierr = mField.set_field_order(0,MBVERTEX,"H1FIELD",1); CHKERRQ(ierr);
  ierr = mField.set_field_order(0,MBTET,"H1FIELD_L2",3); CHKERRQ(ierr);

  FieldCore core2(moab);
  FieldInterface& mField2 = core2;

  //build fields
  ierr = mField2.build_fields(); CHKERRQ(ierr);

  //build finite elements
  ierr = mField2.build_finite_elements(); CHKERRQ(ierr);

  //build adjacencies
  ierr = mField2.build_adjacencies(bit_level0|bit_level1); CHKERRQ(ierr);

  //build problem
  ierr = mField2.build_problems(); CHKERRQ(ierr);

  //partition
  ierr = mField2.partition_problem("PROBLEM_APPROXIMATION"); CHKERRQ(ierr);
  ierr = mField2.partition_finite_elements("PROBLEM_APPROXIMATION"); CHKERRQ(ierr);
  ierr = mField2.partition_ghost_dofs("PROBLEM_APPROXIMATION"); CHKERRQ(ierr);
  ierr = mField2.partition_problem("PROBLEM_APPROXIMATION_REF"); CHKERRQ(ierr);
  ierr = mField2.partition_finite_elements("PROBLEM_APPROXIMATION_REF"); CHKERRQ(ierr);
  ierr = mField2.partition_ghost_dofs("PROBLEM_APPROXIMATION_REF"); CHKERRQ(ierr);

  Vec Dofs_row;
  ierr = mField2.VecCreateGhost("PROBLEM_APPROXIMATION",Row,&Dofs_row); CHKERRQ(ierr);
  Mat Aij;
  ierr = mField2.MatCreateMPIAIJWithArrays("PROBLEM_APPROXIMATION",&Aij); CHKERRQ(ierr);

  struct MyFEMethod: public FEMethod_UpLevelStudent {

    vector<double> g_NTET;

    Mat &Aij;
    Vec& rows_vec;
    MyFEMethod(Interface& _moab,Mat &_Aij,Vec& _rows_vec): FEMethod_UpLevelStudent(_moab,1),Aij(_Aij),rows_vec(_rows_vec) { 
      pcomm = ParallelComm::get_pcomm(&moab,MYPCOMM_INDEX);
    }; 

    ParallelComm* pcomm;
    PetscLogDouble t1,t2;
    PetscLogDouble v1,v2;

    vector<double> g_NTRI;

    double fun(double *fun_coords) { return exp(-fabs(fun_coords[0])/0.05); } //cblas_ddot(3,fun_coords,1,fun_coords,1); }

    PetscErrorCode preProcess() {
      PetscFunctionBegin;
      FEMethod_LowLevelStudent::preProcess();
      g_NTET.resize(4*45);
      ShapeMBTET(&g_NTET[0],G_TET_X45,G_TET_Y45,G_TET_Z45,45);
      PetscSynchronizedPrintf(PETSC_COMM_WORLD,"Start Assemble\n");
      ierr = VecZeroEntries(rows_vec); CHKERRQ(ierr);
      ierr = VecGhostUpdateBegin(rows_vec,INSERT_VALUES,SCATTER_FORWARD); CHKERRQ(ierr);
      ierr = VecGhostUpdateEnd(rows_vec,INSERT_VALUES,SCATTER_FORWARD); CHKERRQ(ierr);
      ierr = MatZeroEntries(Aij); CHKERRQ(ierr);
      ierr = PetscGetTime(&v1); CHKERRQ(ierr);
      ierr = PetscGetCPUTime(&t1); CHKERRQ(ierr);
      PetscFunctionReturn(0);
    }

    PetscErrorCode operator()() {
      PetscFunctionBegin;
      ierr = OpStudentStart_TET(g_NTET); CHKERRQ(ierr);

      g_NTRI.resize(3*4);
      ShapeMBTRI(&g_NTRI[0],G_TRI_X4,G_TRI_Y4,4); 

      SideNumber_multiIndex& side_table = const_cast<SideNumber_multiIndex&>(fe_ent_ptr->get_side_number_table());
      SideNumber_multiIndex::nth_index<1>::type::iterator siit = side_table.get<1>().lower_bound(boost::make_tuple(MBTRI,0));
      SideNumber_multiIndex::nth_index<1>::type::iterator hi_siit = side_table.get<1>().upper_bound(boost::make_tuple(MBTRI,4));
      for(;siit!=hi_siit;siit++) {

	ierr = ShapeFunctions_TRI(siit->ent,g_NTRI);  CHKERRQ(ierr);
	vector< ublas::matrix<FieldData> > FaceNMatrix_nodes;
	ierr = GetGaussRowFaceNMatrix(siit->ent,"H1FIELD",FaceNMatrix_nodes,MBVERTEX); CHKERRQ(ierr);
	//copy(FaceNMatrix_nodes.begin(),FaceNMatrix_nodes.end(),ostream_iterator<ublas::matrix<FieldData> >(cerr," \n")); cerr << endl;

	vector< ublas::matrix<FieldData> > FaceNMatrix_face;
	ierr = GetGaussRowFaceNMatrix(siit->ent,"H1FIELD",FaceNMatrix_face,MBTRI); CHKERRQ(ierr);
	//copy(FaceNMatrix_face.begin(),FaceNMatrix_face.end(),ostream_iterator<ublas::matrix<FieldData> >(cerr," \n")); cerr << endl;

	SideNumber_multiIndex::nth_index<1>::type::iterator siiit = side_table.get<1>().lower_bound(boost::make_tuple(MBEDGE,0));
	SideNumber_multiIndex::nth_index<1>::type::iterator hi_siiit = side_table.get<1>().upper_bound(boost::make_tuple(MBEDGE,6));
	for(;siiit!=hi_siiit;siiit++) {
	  ierr = GetGaussRowFaceNMatrix(siit->ent,"H1FIELD",FaceNMatrix_face,MBEDGE,siiit->ent); CHKERRQ(ierr);
	  //cerr << "ee ";
	  //copy(FaceNMatrix_face.begin(),FaceNMatrix_face.end(),ostream_iterator<ublas::matrix<FieldData> >(cerr," \n")); cerr << endl;
	}


      }
     
      int row_mat = 0;
      vector<vector<DofIdx> > RowGlob(1+6+4+1);
      vector<vector<ublas::matrix<FieldData> > > rowNMatrices(1+6+4+1);
      ierr = GetRowGlobalIndices("H1FIELD",RowGlob[row_mat]); CHKERRQ(ierr);
      ierr = GetGaussRowNMatrix("H1FIELD",rowNMatrices[row_mat++]); CHKERRQ(ierr);
      for(int ee = 0;ee<6;ee++) {
	ierr = GetRowGlobalIndices("H1FIELD",MBEDGE,RowGlob[row_mat],ee); CHKERRQ(ierr);
	if(RowGlob[row_mat].size()!=0) {
	  ierr = GetGaussRowNMatrix("H1FIELD",MBEDGE,rowNMatrices[row_mat++],ee); CHKERRQ(ierr);
	}
      }
      for(int ff = 0;ff<4;ff++) {
	ierr = GetRowGlobalIndices("H1FIELD",MBTRI,RowGlob[row_mat],ff); CHKERRQ(ierr);
	if(RowGlob[row_mat].size()!=0) {
	  ierr = GetGaussRowNMatrix("H1FIELD",MBTRI,rowNMatrices[row_mat++],ff); CHKERRQ(ierr);
	}
      }
      ierr = GetRowGlobalIndices("H1FIELD",MBTET,RowGlob[row_mat]); CHKERRQ(ierr);
      if(RowGlob[row_mat].size() != 0) {
	ierr = GetGaussRowNMatrix("H1FIELD",MBTET,rowNMatrices[row_mat++]); CHKERRQ(ierr);
      }
      int col_mat = 0;
      vector<vector<DofIdx> > ColGlob(1+6+4+1);
      vector<vector<ublas::matrix<FieldData> > > colNMatrices(1+6+4+1);
      ierr = GetColGlobalIndices("H1FIELD",ColGlob[col_mat]); CHKERRQ(ierr);
      ierr = GetGaussColNMatrix("H1FIELD",colNMatrices[col_mat++]); CHKERRQ(ierr);
      for(int ee = 0;ee<6;ee++) {
	ierr = GetColGlobalIndices("H1FIELD",MBEDGE,ColGlob[col_mat],ee); CHKERRQ(ierr);
	if(ColGlob[col_mat].size()!=0) {
	  ierr = GetGaussColNMatrix("H1FIELD",MBEDGE,colNMatrices[col_mat++],ee); CHKERRQ(ierr);
	}
      }
      for(int ff = 0;ff<4;ff++) {
	ierr = GetColGlobalIndices("H1FIELD",MBTRI,ColGlob[col_mat],ff); CHKERRQ(ierr);
	if(ColGlob[col_mat].size()!=0) {
	  ierr = GetGaussColNMatrix("H1FIELD",MBTRI,colNMatrices[col_mat++],ff); CHKERRQ(ierr);
	}
      }
      ierr = GetColGlobalIndices("H1FIELD",MBTET,ColGlob[col_mat]); CHKERRQ(ierr);
      if(ColGlob[col_mat].size() != 0) {
	ierr = GetGaussColNMatrix("H1FIELD",MBTET,colNMatrices[col_mat++]); CHKERRQ(ierr);
      }

      ublas::matrix<FieldData> NTN[row_mat][col_mat];
      ublas::vector<FieldData> F[row_mat];
      int g_dim = g_NTET.size()/4;
      for(int rr = 0;rr<row_mat;rr++) {
	for(int cc = 0;cc<col_mat;cc++) {
	  for(int gg = 0;gg<g_dim;gg++) {
	    ublas::matrix<FieldData> &row_Mat = (rowNMatrices[rr])[gg];
	    ublas::matrix<FieldData> &col_Mat = (colNMatrices[cc])[gg];
	    if(gg == 0) {
	      NTN[rr][cc] = ublas::zero_matrix<FieldData>(row_Mat.size2(),col_Mat.size2());
	    }
	    double w = G_TET_W45[gg];
	    assert(w == w);
	    NTN[rr][cc] += w*V*prod(trans(row_Mat), col_Mat);
	  }
	}
	for(int gg = 0;gg<g_dim;gg++) {
	  ublas::matrix<FieldData> &row_Mat = (rowNMatrices[rr])[gg];
	  if(gg == 0) {
	    F[rr] = ublas::zero_vector<FieldData>(row_Mat.size2());
	  }
	  double val  = fun(&coords_at_Gauss_nodes[gg].data()[0]);
	  F[rr] += V*G_TET_W45[gg]*val*ublas::matrix_row<ublas::matrix<FieldData> >(row_Mat,0);
	}
	if(RowGlob[rr].size()!=F[rr].size()) SETERRQ(PETSC_COMM_SELF,1,"data inconsistency");
	if(RowGlob[rr].size()==0) continue;
	ierr = VecSetValues(rows_vec,RowGlob[rr].size(),&(RowGlob[rr])[0],&(F[rr].data())[0],ADD_VALUES); CHKERRQ(ierr);
	for(int cc = 0;cc<col_mat;cc++) {
	  if(RowGlob[rr].size()!=NTN[rr][cc].size1()) SETERRQ(PETSC_COMM_SELF,1,"data inconsistency");
	  if(ColGlob[cc].size()!=NTN[rr][cc].size2()) SETERRQ(PETSC_COMM_SELF,1,"data inconsistency");
	  if(ColGlob[cc].size()==0) continue;
	  ierr = MatSetValues(Aij,RowGlob[rr].size(),&(RowGlob[rr])[0],ColGlob[cc].size(),&(ColGlob[cc])[0],&(NTN[rr][cc].data())[0],ADD_VALUES); CHKERRQ(ierr);
	}
      }

      ierr = OpStudentEnd(); CHKERRQ(ierr);
      PetscFunctionReturn(0);
    }
    PetscErrorCode postProcess() {
      PetscFunctionBegin;
      ierr = VecGhostUpdateBegin(rows_vec,ADD_VALUES,SCATTER_REVERSE); CHKERRQ(ierr);
      ierr = VecGhostUpdateEnd(rows_vec,ADD_VALUES,SCATTER_REVERSE); CHKERRQ(ierr);
      ierr = VecAssemblyBegin(rows_vec); CHKERRQ(ierr);
      ierr = VecAssemblyEnd(rows_vec); CHKERRQ(ierr);
      ierr = MatAssemblyBegin(Aij,MAT_FINAL_ASSEMBLY); CHKERRQ(ierr);
      ierr = MatAssemblyEnd(Aij,MAT_FINAL_ASSEMBLY); CHKERRQ(ierr);
      ierr = PetscGetTime(&v2);CHKERRQ(ierr);
      ierr = PetscGetCPUTime(&t2);CHKERRQ(ierr);
      PetscSynchronizedPrintf(PETSC_COMM_WORLD,"End Assemble: Rank %d Time = %f CPU Time = %f\n",pcomm->rank(),v2-v1,t2-t1);
      PetscFunctionReturn(0);
    }
  };

  ierr = PetscGetTime(&v2);CHKERRQ(ierr);
  ierr = PetscGetCPUTime(&t2);CHKERRQ(ierr);

  PetscSynchronizedPrintf(PETSC_COMM_WORLD,"Total Rank %d Time = %f CPU Time = %f\n",pcomm->rank(),v2-v1,t2-t1);
  PetscSynchronizedFlush(PETSC_COMM_WORLD);

  MyFEMethod fe_method(moab,Aij,Dofs_row);
  ierr = mField2.loop_finite_elements("PROBLEM_APPROXIMATION","FEAPPROX",fe_method);  CHKERRQ(ierr);
  PetscSynchronizedFlush(PETSC_COMM_WORLD);

  //MatView(Aij,PETSC_VIEWER_DRAW_WORLD);//PETSC_VIEWER_STDOUT_WORLD);
  //std::string wait;
  //std::cin >> wait;

  KSP solver;
  ierr = KSPCreate(PETSC_COMM_WORLD,&solver); CHKERRQ(ierr);
  ierr = KSPSetOperators(solver,Aij,Aij,SAME_NONZERO_PATTERN); CHKERRQ(ierr);
  ierr = KSPSetFromOptions(solver); CHKERRQ(ierr);
  ierr = KSPSetUp(solver); CHKERRQ(ierr);

  Vec solution;
  ierr = VecDuplicate(Dofs_row,&solution); CHKERRQ(ierr);
  ierr = KSPSolve(solver,Dofs_row,solution); CHKERRQ(ierr);
  ierr = VecGhostUpdateBegin(solution,INSERT_VALUES,SCATTER_FORWARD); CHKERRQ(ierr);
  ierr = VecGhostUpdateEnd(solution,INSERT_VALUES,SCATTER_FORWARD); CHKERRQ(ierr);

  ierr = mField2.set_global_VecCreateGhost("PROBLEM_APPROXIMATION",Col,solution,INSERT_VALUES,SCATTER_REVERSE); CHKERRQ(ierr);

  /*struct MyEntMethod: public FieldInterface::EntMethod {
    ErrorCode rval;
    PetscErrorCode ierr;
    Interface& moab;
    
    vector<double> g_NTET;


    Tag th_val;
    MyEntMethod(Interface& _moab): EntMethod(),moab(_moab) {
      double def_VAL = 0;
      rval = moab.tag_get_handle("H1FIELD_VAL",1,MB_TYPE_DOUBLE,th_val,MB_TAG_CREAT|MB_TAG_SPARSE,&def_VAL); CHKERR(rval);
    }

    PetscErrorCode operator()() {
      PetscFunctionBegin;
      if(dof_ptr->get_ent_type()!=MBVERTEX) PetscFunctionReturn(0);
      EntityHandle ent = dof_ptr->get_ent();
      double fval = dof_ptr->get_FieldData();
      rval = moab.tag_set_data(th_val,&ent,1,&fval);  CHKERR_PETSC(rval);
      PetscFunctionReturn(0);
    }

  };
  
  MyEntMethod ent_method(moab);
  ierr = mField2.loop_dofs("PROBLEM_APPROXIMATION","H1FIELD",Col,ent_method); CHKERRQ(ierr);*/

  {
    Tag th_val;
    double def_VAL = 0;
    rval = moab.tag_get_handle("H1FIELD_VAL",1,MB_TYPE_DOUBLE,th_val,MB_TAG_CREAT|MB_TAG_SPARSE,&def_VAL); CHKERR(rval);
    for(_IT_GET_DOFS_FIELD_BY_NAME_FOR_LOOP_(mField,"PROBLEM_APPROXIMATION",dof_ptr)) {
      if(dof_ptr->get_ent_type()!=MBVERTEX) continue;
      EntityHandle ent = dof_ptr->get_ent();
      double fval = dof_ptr->get_FieldData();
      rval = moab.tag_set_data(th_val,&ent,1,&fval);  CHKERR_PETSC(rval);
    }
  }


  ierr = VecDestroy(&solution); CHKERRQ(ierr);
  ierr = KSPDestroy(&solver); CHKERRQ(ierr);
  ierr = VecDestroy(&Dofs_row); CHKERRQ(ierr);
  ierr = MatDestroy(&Aij); CHKERRQ(ierr);

  EntityHandle out_meshset;
  rval = moab.create_meshset(MESHSET_SET,out_meshset); CHKERR_PETSC(rval);
  ierr = mField2.problem_get_FE("PROBLEM_APPROXIMATION","FEAPPROX",out_meshset); CHKERRQ(ierr);
  
  if(pcomm->rank()==0) {
    rval = moab.write_file(mesh_out_file_name,"VTK","",&out_meshset,1); CHKERR_PETSC(rval);
  }

  struct MyFEMethodProject: public FEMethod_UpLevelStudent {
    MyFEMethodProject(Interface& _moab): FEMethod_UpLevelStudent(_moab,1) { 
      pcomm = ParallelComm::get_pcomm(&moab,MYPCOMM_INDEX);
    }; 

    ParallelComm* pcomm;
    PetscLogDouble t1,t2;
    PetscLogDouble v1,v2;

    vector<double> g_NTET;

    PetscErrorCode preProcess() {
      PetscFunctionBegin;
      FEMethod_LowLevelStudent::preProcess();
      g_NTET.resize(4*45);
      ShapeMBTET(&g_NTET[0],G_TET_X45,G_TET_Y45,G_TET_Z45,45);
      PetscSynchronizedPrintf(PETSC_COMM_WORLD,"Start Project\n",pcomm->rank(),v2-v1,t2-t1);
      ierr = PetscGetTime(&v1); CHKERRQ(ierr);
      ierr = PetscGetCPUTime(&t1); CHKERRQ(ierr);
      PetscFunctionReturn(0);
    }
    PetscErrorCode operator()() {
      PetscFunctionBegin;
      ierr = OpStudentStart_TET(g_NTET); CHKERRQ(ierr);
      ierr = ParentData("FEAPPROX"); CHKERRQ(ierr);
      if(ParentMethod==NULL) SETERRQ(PETSC_COMM_SELF,1,"no parent element");
      if(ParentMethod->fe_ent_ptr==NULL) SETERRQ(PETSC_COMM_SELF,1,"no parent element");
      //
      EntityHandle parent_ent = ParentMethod->fe_ent_ptr->get_ent();
      int num_nodes;
      rval = moab.get_connectivity(parent_ent,ParentMethod->conn,num_nodes,true); CHKERR_PETSC(rval);
      assert(num_nodes == 4);
      ParentMethod->coords.resize(12);
      rval = moab.get_coords(ParentMethod->conn,num_nodes,&ParentMethod->coords[0]); CHKERR_PETSC(rval);
      int g_dim = g_NTET.size()/4;
      vector<double> G_X(g_dim),G_Y(g_dim),G_Z(g_dim);
      ublas::vector<double,ublas::bounded_array<double, 3> > glob_coords;
      ublas::vector<double,ublas::bounded_array<double, 3> > loc_coords;
      double X[1] = { 0 },Y[1] = { 0 },Z[1] = { 0 };
      double NTET0[4];
      ShapeMBTET(NTET0,X,Y,Z,1);
      for(int gg = 0;gg<g_dim;gg++) {
	glob_coords.resize(3);
	loc_coords.resize(3);
	glob_coords[0] = cblas_ddot(4,&g_NTET[gg*4],1,&coords[0],3);
	glob_coords[1] = cblas_ddot(4,&g_NTET[gg*4],1,&coords[1],3);
	glob_coords[2] = cblas_ddot(4,&g_NTET[gg*4],1,&coords[2],3);
	ierr = ShapeMBTET_inverse(NTET0,diffNTET,&ParentMethod->coords[0],&glob_coords.data()[0],&loc_coords.data()[0]); CHKERRQ(ierr);
	assert(loc_coords[0]>=0);
	assert(loc_coords[0]<=1);
	G_X[gg] = loc_coords[0];
	assert(loc_coords[1]>=0);
	assert(loc_coords[1]<=1);
	G_Y[gg] = loc_coords[1];
	assert(loc_coords[2]>=0);
	assert(loc_coords[2]<=1);
	G_Z[gg] = loc_coords[2];
      }
      vector<double> ParentMethod_g_NTET(4*g_dim);
      ShapeMBTET(&ParentMethod_g_NTET[0],&G_X[0],&G_Y[0],&G_Z[0],g_dim);
      ierr = ParentMethod->InitDataStructures(); CHKERRQ(ierr);
      ierr = ParentMethod->DataOp(); CHKERRQ(ierr);
      ierr = ParentMethod->ShapeFunctions_TET(ParentMethod_g_NTET);
      ierr = ParentMethod->Data_at_GaussPoints(); CHKERRQ(ierr);
      Data_at_Gauss_pt &parent_data_at_gauss_pt = ParentMethod->data_at_gauss_pt;
      Data_at_Gauss_pt::iterator diit = parent_data_at_gauss_pt.find("H1FIELD");
      if(diit==parent_data_at_gauss_pt.end()) SETERRQ(PETSC_COMM_SELF,1,"no H1FIELD !!!");
      vector<ublas::vector<FieldData> > &data = diit->second;
      //
      vector<ublas::matrix<FieldData> > RowN;
      vector<ublas::matrix<FieldData> > ColN;
      ierr = GetGaussRowNMatrix("H1FIELD_L2",MBTET,RowN); CHKERRQ(ierr);
      ierr = GetGaussColNMatrix("H1FIELD_L2",MBTET,ColN); CHKERRQ(ierr);
      ublas::matrix<FieldData> NTN;
      NTN = ublas::zero_matrix<FieldData>(RowN[0].size2(),ColN[0].size2());
      ublas::vector<FieldData> F;
      F = ublas::zero_vector<FieldData>(RowN[0].size2());
      for(int gg = 0;gg<g_dim;gg++) {
	double w = G_TET_W45[gg];
	assert(w == w);
	NTN += w*V*prod(trans(RowN[gg]), ColN[gg]);
	double val  = data[gg][0];
	F += V*w*val*ublas::matrix_row<ublas::matrix<FieldData> >(RowN[gg],0);
      }
      ublas::matrix<FieldData> L(NTN.size1(),NTN.size2());
      ublas::vector<FieldData> x = F;
      cholesky_decompose(NTN,L);
      cholesky_solve(L,x,ublas::lower());
      FENumeredDofMoFEMEntity_multiIndex::index<Unique_mi_tag>::type::iterator 
	ddiit = row_multiIndex->get<Unique_mi_tag>().begin(),hi_ddiit = row_multiIndex->get<Unique_mi_tag>().end();
      unsigned int dd = 0;
      for(;ddiit!=hi_ddiit;ddiit++,dd++) {
	FieldData &val = const_cast<FieldData&>(ddiit->get_FieldData());
	val = x[ddiit->get_EntDofIdx()];
      }
      if(dd!=F.size()) SETERRQ(PETSC_COMM_SELF,1,"data inconsistency");
      ierr = OpStudentEnd(); CHKERRQ(ierr);
      PetscFunctionReturn(0);
    }
    PetscErrorCode postProcess() {
      PetscFunctionBegin;
      ierr = PetscGetTime(&v2); CHKERRQ(ierr);
      ierr = PetscGetCPUTime(&t2); CHKERRQ(ierr);
      PetscSynchronizedPrintf(PETSC_COMM_WORLD,"End Project: Rank %d Time = %f CPU Time = %f\n",pcomm->rank(),v2-v1,t2-t1);
      PetscFunctionReturn(0);
    }

  };
  
  MyFEMethodProject fe_proj_method(moab);
  ierr = mField2.loop_finite_elements("PROBLEM_APPROXIMATION_REF","FEAPPROX_L2",fe_proj_method);  CHKERRQ(ierr);
  PetscSynchronizedFlush(PETSC_COMM_WORLD);

  struct MyFEMethodPostProc: public FEMethod_UpLevelStudent {
    vector<double> g_NTET;


    Tag th_val;
    MyFEMethodPostProc(Interface& _moab): FEMethod_UpLevelStudent(_moab,1),moab_post_proc(mb_instance_post_proc) { 
      pcomm = ParallelComm::get_pcomm(&moab,MYPCOMM_INDEX);
      double def_VAL = 0;
      rval = moab_post_proc.tag_get_handle("H1FIELD_L2_VAL",1,MB_TYPE_DOUBLE,th_val,MB_TAG_CREAT|MB_TAG_SPARSE,&def_VAL); CHKERR_THROW(rval);
    }; 

    Interface& moab_post_proc;
    Core mb_instance_post_proc;

    ParallelComm* pcomm;
    PetscLogDouble t1,t2;
    PetscLogDouble v1,v2;

    PetscErrorCode preProcess() {
      PetscFunctionBegin;
      g_NTET.resize(4*10);
      double GX[10] = { 0,1,0,0, 0.5*(GX[0]+GX[1]), 0.5*(GX[1]+GX[2]), 0.5*(GX[2]+GX[0]), 0.5*(GX[0]+GX[3]), 0.5*(GX[1]+GX[3]), 0.5*(GX[2]+GX[3]) };
      double GY[10] = { 0,0,1,0, 0.5*(GY[0]+GY[1]), 0.5*(GY[1]+GY[2]), 0.5*(GY[2]+GY[0]), 0.5*(GY[0]+GY[3]), 0.5*(GY[1]+GY[3]), 0.5*(GY[2]+GY[3]) };
      double GZ[10] = { 0,0,0,1, 0.5*(GZ[0]+GZ[1]), 0.5*(GZ[1]+GZ[2]), 0.5*(GZ[2]+GZ[0]), 0.5*(GZ[0]+GZ[3]), 0.5*(GZ[1]+GZ[3]), 0.5*(GZ[2]+GZ[3]) };
      ShapeMBTET(&g_NTET[0],GX,GY,GZ,10);
      PetscSynchronizedPrintf(PETSC_COMM_WORLD,"Start PostProc\n",pcomm->rank(),v2-v1,t2-t1);
      ierr = PetscGetTime(&v1); CHKERRQ(ierr);
      ierr = PetscGetCPUTime(&t1); CHKERRQ(ierr);
      FEMethod_LowLevelStudent::preProcess();
      PetscFunctionReturn(0);
    }
    PetscErrorCode operator()() {
      PetscFunctionBegin;
      ierr = OpStudentStart_TET(g_NTET); CHKERRQ(ierr);

      Data_at_Gauss_pt::iterator diit = data_at_gauss_pt.find("H1FIELD_L2");
      if(diit==data_at_gauss_pt.end()) SETERRQ(PETSC_COMM_SELF,1,"no H1FIELD L2 !!!");
      vector< ublas::vector<FieldData> > &data = diit->second;
      if(data.size()!=10) SETERRQ(PETSC_COMM_SELF,1,"data inconsistency");

      EntityHandle nodes[4+6];
      for(int nn = 0;nn<4;nn++) {
	rval = moab_post_proc.create_vertex(&coords[3*nn],nodes[nn]); CHKERR_PETSC(rval);
      }
      double coord[3];
      fill(&coord[0],&coord[3],0);
      cblas_daxpy(3,0.5,&coords[3*0],1,coord,1);
      cblas_daxpy(3,0.5,&coords[3*1],1,coord,1);
      rval = moab_post_proc.create_vertex(coord,nodes[4]); CHKERR_PETSC(rval);
      fill(&coord[0],&coord[3],0);
      cblas_daxpy(3,0.5,&coords[3*1],1,coord,1);
      cblas_daxpy(3,0.5,&coords[3*2],1,coord,1);
      rval = moab_post_proc.create_vertex(coord,nodes[5]); CHKERR_PETSC(rval);
      fill(&coord[0],&coord[3],0);
      cblas_daxpy(3,0.5,&coords[3*2],1,coord,1);
      cblas_daxpy(3,0.5,&coords[3*0],1,coord,1);
      rval = moab_post_proc.create_vertex(coord,nodes[6]); CHKERR_PETSC(rval);
      fill(&coord[0],&coord[3],0);
      cblas_daxpy(3,0.5,&coords[3*0],1,coord,1);
      cblas_daxpy(3,0.5,&coords[3*3],1,coord,1);
      rval = moab_post_proc.create_vertex(coord,nodes[7]); CHKERR_PETSC(rval);
      fill(&coord[0],&coord[3],0);
      cblas_daxpy(3,0.5,&coords[3*1],1,coord,1);
      cblas_daxpy(3,0.5,&coords[3*3],1,coord,1);
      rval = moab_post_proc.create_vertex(coord,nodes[8]); CHKERR_PETSC(rval);
      fill(&coord[0],&coord[3],0);
      cblas_daxpy(3,0.5,&coords[3*2],1,coord,1);
      cblas_daxpy(3,0.5,&coords[3*3],1,coord,1);
      rval = moab_post_proc.create_vertex(coord,nodes[9]); CHKERR_PETSC(rval);
      EntityHandle tets[4*8];
      tet_type_6(moab_post_proc,nodes,&nodes[4],tets);
      EntityHandle ref_tets[8];
      for(int tt = 0;tt<8;tt++) {
	rval = moab_post_proc.create_element(MBTET,&tets[4*tt],4,ref_tets[tt]); CHKERR_PETSC(rval);
      }
      for(int nn = 0;nn<10;nn++) {
	double vval = (data[nn])[0];
	rval = moab_post_proc.tag_set_data(th_val,&nodes[nn],1,&vval); CHKERR_PETSC(rval);
      }

      ierr = OpStudentEnd(); CHKERRQ(ierr);
      PetscFunctionReturn(0);
    }
    PetscErrorCode postProcess() {
      PetscFunctionBegin;
      ierr = PetscGetTime(&v2); CHKERRQ(ierr);
      ierr = PetscGetCPUTime(&t2); CHKERRQ(ierr);
      PetscSynchronizedPrintf(PETSC_COMM_WORLD,"End PostProc: Rank %d Time = %f CPU Time = %f\n",pcomm->rank(),v2-v1,t2-t1);
      ParallelComm* pcomm_post_proc = ParallelComm::get_pcomm(&moab_post_proc,MYPCOMM_INDEX);
      if(pcomm_post_proc == NULL) pcomm_post_proc =  new ParallelComm(&moab_post_proc,PETSC_COMM_WORLD);
      for(unsigned int rr = 1; rr<pcomm_post_proc->size();rr++) {
	Range tets;
	rval = moab_post_proc.get_entities_by_type(0,MBTET,tets); CHKERR_PETSC(rval);
	rval = pcomm_post_proc->broadcast_entities(rr,tets); CHKERR(rval);

      }
      PetscFunctionReturn(0);
    }

  };

  MyFEMethodPostProc fe_post_proc_method(moab);
  ierr = mField2.loop_finite_elements("PROBLEM_APPROXIMATION_REF","FEAPPROX_L2",fe_post_proc_method);  CHKERRQ(ierr);
  PetscSynchronizedFlush(PETSC_COMM_WORLD);
  if(pcomm->rank()==0) {
    rval = fe_post_proc_method.moab_post_proc.write_file("out_post_proc.vtk","VTK",""); CHKERR_PETSC(rval);
  }

  PetscFinalize();

}

