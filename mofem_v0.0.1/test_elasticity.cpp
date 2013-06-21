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
#include "moabFEMethod_Student.hpp"
#include "cholesky.hpp"
#include <petscksp.h>

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

  moabField_Core core(moab);
  moabField& mField = core;

  Range CubitSideSets_meshsets;
  ierr = mField.get_CubitBCType_meshsets(SideSet,CubitSideSets_meshsets); CHKERRQ(ierr);

  //ref meshset ref level 0
  ierr = mField.seed_ref_level_3D(0,0); CHKERRQ(ierr);
  BitRefLevel bit_level0;
  bit_level0.set(0);
  EntityHandle meshset_level0;
  rval = moab.create_meshset(MESHSET_SET,meshset_level0); CHKERR_PETSC(rval);
  ierr = mField.seed_ref_level_3D(0,bit_level0); CHKERRQ(ierr);
  ierr = mField.refine_get_ents(bit_level0,meshset_level0); CHKERRQ(ierr);

  //Fields
  ierr = mField.add_BitFieldId("DISPLACEMENT",H1,3); CHKERRQ(ierr);

  //FE
  ierr = mField.add_MoFEMFE("ELASTIC"); CHKERRQ(ierr);

  //Define rows/cols and element data
  ierr = mField.modify_MoFEMFE_row_add_bit("ELASTIC","DISPLACEMENT"); CHKERRQ(ierr);
  ierr = mField.modify_MoFEMFE_col_add_bit("ELASTIC","DISPLACEMENT"); CHKERRQ(ierr);
  ierr = mField.modify_MoFEMFE_data_add_bit("ELASTIC","DISPLACEMENT"); CHKERRQ(ierr);

  //define problems
  ierr = mField.add_BitProblemId("ELASTIC_MECHANICS"); CHKERRQ(ierr);

  //set finite elements for problems
  ierr = mField.modify_problem_MoFEMFE_add_bit("ELASTIC_MECHANICS","ELASTIC"); CHKERRQ(ierr);

  //set refinment level for problem
  ierr = mField.modify_problem_ref_level_add_bit("ELASTIC_MECHANICS",bit_level0); CHKERRQ(ierr);

  //add entitities (by tets) to the field
  ierr = mField.add_ents_to_field_by_TETs(0,"DISPLACEMENT"); CHKERRQ(ierr);

  //add finite elements entities
  ierr = mField.add_ents_to_MoFEMFE_by_bit_ref(bit_level0,"ELASTIC",MBTET); CHKERRQ(ierr);

  //set app. order
  ierr = mField.set_field_order(0,MBTET,"DISPLACEMENT",1); CHKERRQ(ierr);
  ierr = mField.set_field_order(0,MBTRI,"DISPLACEMENT",1); CHKERRQ(ierr);
  ierr = mField.set_field_order(0,MBEDGE,"DISPLACEMENT",1); CHKERRQ(ierr);
  ierr = mField.set_field_order(0,MBVERTEX,"DISPLACEMENT",1); CHKERRQ(ierr);

  //build field
  ierr = mField.build_fields(); CHKERRQ(ierr);

  //build finite elemnts
  ierr = mField.build_finite_elements(); CHKERRQ(ierr);

  //build adjacencies
  ierr = mField.build_adjacencies(bit_level0); CHKERRQ(ierr);

  //build problem
  ierr = mField.build_problems(); CHKERRQ(ierr);

  //partition
  ierr = mField.partition_problems("ELASTIC_MECHANICS"); CHKERRQ(ierr);
  ierr = mField.partition_finite_elements("ELASTIC_MECHANICS"); CHKERRQ(ierr);
  ierr = mField.partition_ghost_dofs("ELASTIC_MECHANICS"); CHKERRQ(ierr);

  //create matrices
  Vec F;
  ierr = mField.VecCreateGhost("ELASTIC_MECHANICS",Row,&F); CHKERRQ(ierr);
  Mat Aij;
  ierr = mField.MatCreateMPIAIJWithArrays("ELASTIC_MECHANICS",&Aij); CHKERRQ(ierr);

  struct ElasticFEMethod: public FEMethod_Student {

    Mat &Aij;
    Vec& F;
    ElasticFEMethod(
      Interface& _moab,Mat &_Aij,Vec& _F,
      double _lambda,double _mu,Range &_SideSet1,Range &_SideSet2): 
      FEMethod_Student(_moab,1),Aij(_Aij),F(_F),
      lambda(_lambda),mu(_mu),
      SideSet1(_SideSet1),SideSet2(_SideSet2) { 
      pcomm = ParallelComm::get_pcomm(&moab,MYPCOMM_INDEX);

      RowGlob.resize(1+6+4+1);
      rowNMatrices.resize(1+6+4+1);
      rowDiffNMatrices.resize(1+6+4+1);
      rowBMatrices.resize(1+6+4+1);
      ColGlob.resize(1+6+4+1);
      colNMatrices.resize(1+6+4+1);
      colDiffNMatrices.resize(1+6+4+1);
      colBMatrices.resize(1+6+4+1);

      //
      Range SideSet1Edges,SideSet1Nodes;
      rval = moab.get_adjacencies(SideSet1,1,false,SideSet1Edges,Interface::UNION); CHKERR_THROW(rval);
      rval = moab.get_connectivity(SideSet1,SideSet1Nodes,true); CHKERR_THROW(rval);
      SideSet1.insert(SideSet1Edges.begin(),SideSet1Edges.end());
      SideSet1.insert(SideSet1Nodes.begin(),SideSet1Nodes.end());

    }; 

    ErrorCode rval;
    
    ParallelComm* pcomm;
    PetscLogDouble t1,t2;
    PetscLogDouble v1,v2;

    double lambda,mu;
    ublas::matrix<FieldData> D_lambda,D_mu,D;

    Range& SideSet1;
    Range& SideSet2;

    int row_mat,col_mat;
    vector<vector<DofIdx> > RowGlob;
    vector<vector<ublas::matrix<FieldData> > > rowNMatrices;
    vector<vector<ublas::matrix<FieldData> > > rowDiffNMatrices;
    vector<vector<ublas::matrix<FieldData> > > rowBMatrices;
    vector<vector<DofIdx> > ColGlob;
    vector<vector<ublas::matrix<FieldData> > > colNMatrices;
    vector<vector<ublas::matrix<FieldData> > > colDiffNMatrices;
    vector<vector<ublas::matrix<FieldData> > > colBMatrices;

    vector<DofIdx> DirihletBC;
    vector<FieldData> DirihletBCDiagVal;
    
    PetscErrorCode preProcess() {
      PetscFunctionBegin;
      FEMethod_Core::preProcess();
      PetscSynchronizedPrintf(PETSC_COMM_WORLD,"Start Assembly\n",pcomm->rank(),v2-v1,t2-t1);
      ierr = PetscGetTime(&v1); CHKERRQ(ierr);
      ierr = PetscGetCPUTime(&t1); CHKERRQ(ierr);
      g_NTET.resize(4*4);
      ShapeMBTET(&g_NTET[0],G_TET_X4,G_TET_Y4,G_TET_Z4,4);
      ierr = VecZeroEntries(F); CHKERRQ(ierr);
      ierr = VecGhostUpdateBegin(F,INSERT_VALUES,SCATTER_FORWARD); CHKERRQ(ierr);
      ierr = VecGhostUpdateEnd(F,INSERT_VALUES,SCATTER_FORWARD); CHKERRQ(ierr);
      ierr = MatZeroEntries(Aij); CHKERRQ(ierr);
      // See FEAP - - A Finite Element Analysis Program
      D_lambda = ublas::zero_matrix<FieldData>(6,6);
      for(int rr = 0;rr<3;rr++) {
	ublas::matrix_row<ublas::matrix<FieldData> > row_D_lambda(D_lambda,rr);
	for(int cc = 0;cc<3;cc++) {
	  row_D_lambda[cc] = 1;
	}
      }
      D_mu = ublas::zero_matrix<FieldData>(6,6);
      for(int rr = 0;rr<3;rr++) {
	D_mu(rr,rr) = rr<3 ? 2 : 1;
      }
      D = lambda*D_lambda + mu*D_mu;

      PetscFunctionReturn(0);
    }
    PetscErrorCode postProcess() {
      PetscFunctionBegin;
      ierr = VecGhostUpdateBegin(F,ADD_VALUES,SCATTER_REVERSE); CHKERRQ(ierr);
      ierr = VecGhostUpdateEnd(F,ADD_VALUES,SCATTER_REVERSE); CHKERRQ(ierr);
      ierr = VecAssemblyBegin(F); CHKERRQ(ierr);
      ierr = VecAssemblyEnd(F); CHKERRQ(ierr);
      ierr = MatAssemblyBegin(Aij,MAT_FINAL_ASSEMBLY); CHKERRQ(ierr);
      ierr = MatAssemblyEnd(Aij,MAT_FINAL_ASSEMBLY); CHKERRQ(ierr);
      ierr = PetscGetTime(&v2); CHKERRQ(ierr);
      ierr = PetscGetCPUTime(&t2); CHKERRQ(ierr);
      PetscSynchronizedPrintf(PETSC_COMM_WORLD,"End Assembly: Rank %d Time = %f CPU Time = %f\n",pcomm->rank(),v2-v1,t2-t1);
      PetscFunctionReturn(0);
    }
    
    PetscErrorCode GetMatrices() {
      PetscFunctionBegin;
      //indicies ROWS
      row_mat = 0;
      ierr = GetRowIndices("DISPLACEMENT",RowGlob[row_mat]); CHKERRQ(ierr);
      ierr = GetGaussRowNMatrix("DISPLACEMENT",rowNMatrices[row_mat]); CHKERRQ(ierr);
      ierr = GetGaussRowDiffNMatrix("DISPLACEMENT",rowDiffNMatrices[row_mat]); CHKERRQ(ierr);
      ierr = MakeBMatrix3D("DISPLACEMENT",rowDiffNMatrices[row_mat],rowBMatrices[row_mat]);  CHKERRQ(ierr);
      row_mat++;
      for(int ee = 0;ee<6;ee++) { //edges matrices
	ierr = GetRowIndices("DISPLACEMENT",MBEDGE,RowGlob[row_mat],ee); CHKERRQ(ierr);
	if(RowGlob[row_mat].size()!=0) {
	  ierr = GetGaussRowNMatrix("DISPLACEMENT",MBEDGE,rowNMatrices[row_mat],ee); CHKERRQ(ierr);
	  ierr = GetGaussRowDiffNMatrix("DISPLACEMENT",MBEDGE,rowDiffNMatrices[row_mat],ee); CHKERRQ(ierr);
	  ierr = MakeBMatrix3D("DISPLACEMENT",rowDiffNMatrices[row_mat],rowBMatrices[row_mat]);  CHKERRQ(ierr);
	  row_mat++;
	}
      }
      for(int ff = 0;ff<4;ff++) { //faces matrices
	ierr = GetRowIndices("DISPLACEMENT",MBTRI,RowGlob[row_mat],ff); CHKERRQ(ierr);
	if(RowGlob[row_mat].size()!=0) {
	  ierr = GetGaussRowNMatrix("DISPLACEMENT",MBTRI,rowNMatrices[row_mat],ff); CHKERRQ(ierr);
	  ierr = GetGaussRowDiffNMatrix("DISPLACEMENT",MBTRI,rowDiffNMatrices[row_mat],ff); CHKERRQ(ierr);
	  ierr = MakeBMatrix3D("DISPLACEMENT",rowDiffNMatrices[row_mat],rowBMatrices[row_mat]);  CHKERRQ(ierr);
	  row_mat++;
	}
      }
      ierr = GetRowIndices("DISPLACEMENT",MBTET,RowGlob[row_mat]); CHKERRQ(ierr);
      if(RowGlob[row_mat].size() != 0) { //volume matrices
	ierr = GetGaussRowNMatrix("DISPLACEMENT",MBTET,rowNMatrices[row_mat++]); CHKERRQ(ierr);
	ierr = GetGaussRowDiffNMatrix("DISPLACEMENT",MBTET,rowDiffNMatrices[row_mat]); CHKERRQ(ierr);
	ierr = MakeBMatrix3D("DISPLACEMENT",rowDiffNMatrices[row_mat],rowBMatrices[row_mat]);  CHKERRQ(ierr);
	row_mat++;
      }

      //indicies COLS
      col_mat = 0;
      ierr = GetColIndices("DISPLACEMENT",ColGlob[col_mat]); CHKERRQ(ierr);
      ierr = GetGaussColNMatrix("DISPLACEMENT",colNMatrices[col_mat]); CHKERRQ(ierr);
      ierr = GetGaussColDiffNMatrix("DISPLACEMENT",colDiffNMatrices[col_mat]); CHKERRQ(ierr);
      ierr = MakeBMatrix3D("DISPLACEMENT",colDiffNMatrices[col_mat],colBMatrices[col_mat]);  CHKERRQ(ierr);
      col_mat++;
      for(int ee = 0;ee<6;ee++) { //edges matrices
	ierr = GetColIndices("DISPLACEMENT",MBEDGE,ColGlob[col_mat],ee); CHKERRQ(ierr);
	if(ColGlob[col_mat].size()!=0) {
	  ierr = GetGaussColNMatrix("DISPLACEMENT",MBEDGE,colNMatrices[col_mat],ee); CHKERRQ(ierr);
	  ierr = GetGaussColDiffNMatrix("DISPLACEMENT",MBEDGE,colDiffNMatrices[col_mat],ee); CHKERRQ(ierr);
	  ierr = MakeBMatrix3D("DISPLACEMENT",colDiffNMatrices[col_mat],colBMatrices[col_mat]);  CHKERRQ(ierr);
	  col_mat++;
	}
      }
      for(int ff = 0;ff<4;ff++) { //faces matrices
	ierr = GetColIndices("DISPLACEMENT",MBTRI,ColGlob[col_mat],ff); CHKERRQ(ierr);
	if(ColGlob[col_mat].size()!=0) {
	  ierr = GetGaussColNMatrix("DISPLACEMENT",MBTRI,colNMatrices[col_mat],ff); CHKERRQ(ierr);
	  ierr = GetGaussColDiffNMatrix("DISPLACEMENT",MBTRI,colDiffNMatrices[col_mat],ff); CHKERRQ(ierr);
	  ierr = MakeBMatrix3D("DISPLACEMENT",colDiffNMatrices[col_mat],colBMatrices[col_mat]);  CHKERRQ(ierr);
	  col_mat++;
	}
      }
      ierr = GetColIndices("DISPLACEMENT",MBTET,ColGlob[col_mat]); CHKERRQ(ierr);
      if(ColGlob[col_mat].size() != 0) { //volume matrices
	ierr = GetGaussColNMatrix("DISPLACEMENT",MBTET,colNMatrices[col_mat++]); CHKERRQ(ierr);
	ierr = GetGaussColDiffNMatrix("DISPLACEMENT",MBTET,colDiffNMatrices[col_mat]); CHKERRQ(ierr);
	ierr = MakeBMatrix3D("DISPLACEMENT",colDiffNMatrices[col_mat],colBMatrices[col_mat]);  CHKERRQ(ierr);
	col_mat++;
      }

      //Boundary Condition
      //Dirihlet form SideSet1
      DirihletBC.resize(0);
      Range::iterator siit1 = SideSet1.begin();
      for(;siit1!=SideSet1.end();siit1++) {
	FENumeredDofMoFEMEntity_multiIndex::index<MoABEnt_mi_tag>::type::iterator riit = row_multiIndex->get<MoABEnt_mi_tag>().lower_bound(*siit1);
	FENumeredDofMoFEMEntity_multiIndex::index<MoABEnt_mi_tag>::type::iterator hi_riit = row_multiIndex->get<MoABEnt_mi_tag>().upper_bound(*siit1);
	for(;riit!=hi_riit;riit++) {
	  if(riit->get_name()!="DISPLACEMENT") continue;
	  // all fixed
	  // if some ranks are selected then we could apply BC in particular direction
	  DirihletBC.push_back(riit->get_petsc_gloabl_dof_idx());
	  for(int cc = 0;cc<col_mat;cc++) {
	    vector<DofIdx>::iterator it = find(ColGlob[cc].begin(),ColGlob[cc].end(),DirihletBC.back());
	    if( it!=ColGlob[cc].end() ) *it = -1; // of idx is set -1 column is not assembled
	  }
	  for(int rr = 0;rr<col_mat;rr++) {
	    vector<DofIdx>::iterator it = find(RowGlob[rr].begin(),RowGlob[rr].end(),DirihletBC.back());
	    if( it!=RowGlob[rr].end() ) *it = -1; // of idx is set -1 row is not assembled
	  }
	}
      }
      DirihletBCDiagVal.resize(DirihletBC.size());
      fill(DirihletBCDiagVal.begin(),DirihletBCDiagVal.end(),1);

      PetscFunctionReturn(0);
    }

    PetscErrorCode operator()() {
      PetscFunctionBegin;
      ierr = OpStudentStart(); CHKERRQ(ierr);
      ierr = GetMatrices(); CHKERRQ(ierr);

      ublas::matrix<FieldData> K[row_mat][col_mat];
      ublas::vector<FieldData> f[row_mat];
      int g_dim = g_NTET.size()/4;
      for(int rr = 0;rr<row_mat;rr++) {
	for(int cc = 0;cc<col_mat;cc++) {
	  for(int gg = 0;gg<g_dim;gg++) {
	    ublas::matrix<FieldData> &row_Mat = (rowBMatrices[rr])[gg];
	    ublas::matrix<FieldData> &col_Mat = (rowBMatrices[cc])[gg];
	    ///K matrices
	    if(gg == 0) {
	      K[rr][cc] = ublas::zero_matrix<FieldData>(row_Mat.size2(),col_Mat.size2());
	    }
	    double w = V*G_TET_W45[gg];
	    ublas::matrix<FieldData> BTD = prod( trans(row_Mat), w*D );
	    K[rr][cc] += prod(BTD , col_Mat ); // int BT*D*B
	  }
	}
	for(int gg = 0;gg<g_dim;gg++) {
	  ublas::matrix<FieldData> &row_Mat = (rowBMatrices[rr])[gg];
	  if(gg == 0) f[rr] = ublas::zero_vector<FieldData>(row_Mat.size2());
	  //ublas::matrix<FieldData> &row_Mat = (rowDiffNMatrices[rr])[gg];
	  ///f matrices

	}
	if(RowGlob[rr].size()!=f[rr].size()) SETERRQ(PETSC_COMM_SELF,1,"data inconsistency");
	if(RowGlob[rr].size()==0) continue;
	ierr = VecSetValues(F,RowGlob[rr].size(),&(RowGlob[rr])[0],&(f[rr].data())[0],ADD_VALUES); CHKERRQ(ierr);
	for(int cc = 0;cc<col_mat;cc++) {
	  if(ColGlob[cc].size()==0) continue;
	  if(RowGlob[rr].size()!=K[rr][cc].size1()) SETERRQ(PETSC_COMM_SELF,1,"data inconsistency");
	  if(ColGlob[cc].size()!=K[rr][cc].size2()) SETERRQ(PETSC_COMM_SELF,1,"data inconsistency");
	  ierr = MatSetValues(Aij,RowGlob[rr].size(),&(RowGlob[rr])[0],ColGlob[cc].size(),&(ColGlob[cc])[0],&(K[rr][cc].data())[0],ADD_VALUES); CHKERRQ(ierr);
	}
      }
      //DirihletBC -> values on diagonal
      ierr = MatSetValues(Aij,DirihletBC.size(),&DirihletBC[0],DirihletBC.size(),&DirihletBC[0],&DirihletBCDiagVal[0],ADD_VALUES); CHKERRQ(ierr);

      ierr = OpStudentEnd(); CHKERRQ(ierr);
      PetscFunctionReturn(0);
    }

  };

  Range SideSet1,SideSet2;
  ierr = mField.get_Cubit_msId_entities_by_dimension(1,SideSet,2,SideSet1,true); CHKERRQ(ierr);
  ierr = mField.get_Cubit_msId_entities_by_dimension(2,SideSet,2,SideSet2,true); CHKERRQ(ierr);
  PetscPrintf(PETSC_COMM_WORLD,"Nb. faces in SideSet 1 : %u\n",SideSet1.size());
  PetscPrintf(PETSC_COMM_WORLD,"Nb. faces in SideSet 1 : %u\n",SideSet2.size());

  ElasticFEMethod MyFE(moab,Aij,F,1,1,SideSet1,SideSet2);
  ierr = mField.loop_finite_elements("ELASTIC_MECHANICS","ELASTIC",MyFE);  CHKERRQ(ierr);
  PetscSynchronizedFlush(PETSC_COMM_WORLD);

  //detroy matrices
  VecDestroy(&F);
  MatDestroy(&Aij);

  ierr = PetscGetTime(&v2);CHKERRQ(ierr);
  ierr = PetscGetCPUTime(&t2);CHKERRQ(ierr);

  PetscSynchronizedPrintf(PETSC_COMM_WORLD,"Total Rank %d Time = %f CPU Time = %f\n",pcomm->rank(),v2-v1,t2-t1);
  PetscSynchronizedFlush(PETSC_COMM_WORLD);

  PetscFinalize();

}

