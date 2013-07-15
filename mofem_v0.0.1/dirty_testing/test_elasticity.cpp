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
#include "cholesky.hpp"
#include <petscksp.h>

/* 
 * From: Demkowicz
 * To discuss basic a priori error estimates, we return now to basic mathemat-
 * ical issues. We have learned from Cea’s lemma that, for the coercive case,
 * the actual approximation error can always be bounded by the product of a
 * mesh-independent constant and the best approximation error. Thus, to es-
 * timate the approximation error, it is sufficient to estimate the best
 * approxi- mation error. By definition, the best approximation error is always
 * bounded by the norm of the difference between the exact solution and any
 * particular choice of a function that “lives” in the FE space. */


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

  PetscBool flg = PETSC_TRUE;
  char mesh_file_name[255];
  ierr = PetscOptionsGetString(PETSC_NULL,"-my_file",mesh_file_name,255,&flg); CHKERRQ(ierr);
  if(flg != PETSC_TRUE) {
    SETERRQ(PETSC_COMM_SELF,1,"*** ERROR -my_file (MESH FILE NEEDED)");
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
  ierr = mField.add_field("DISPLACEMENT",H1,3); CHKERRQ(ierr);

  //FE
  ierr = mField.add_finite_element("ELASTIC"); CHKERRQ(ierr);

  //Define rows/cols and element data
  ierr = mField.modify_finite_element_add_field_row("ELASTIC","DISPLACEMENT"); CHKERRQ(ierr);
  ierr = mField.modify_finite_element_add_field_col("ELASTIC","DISPLACEMENT"); CHKERRQ(ierr);
  ierr = mField.modify_finite_element_add_field_data("ELASTIC","DISPLACEMENT"); CHKERRQ(ierr);

  //define problems
  ierr = mField.add_problem("ELASTIC_MECHANICS"); CHKERRQ(ierr);

  //set finite elements for problems
  ierr = mField.modify_problem_add_finite_element("ELASTIC_MECHANICS","ELASTIC"); CHKERRQ(ierr);

  //set refinment level for problem
  ierr = mField.modify_problem_ref_level_add_bit("ELASTIC_MECHANICS",bit_level0); CHKERRQ(ierr);

  //add entitities (by tets) to the field
  ierr = mField.add_ents_to_field_by_TETs(0,"DISPLACEMENT"); CHKERRQ(ierr);

  //add finite elements entities
  ierr = mField.add_ents_to_finite_element_EntType_by_bit_ref(bit_level0,"ELASTIC",MBTET); CHKERRQ(ierr);

  //set app. order
  ierr = mField.set_field_order(0,MBTET,"DISPLACEMENT",4); CHKERRQ(ierr);
  ierr = mField.set_field_order(0,MBTRI,"DISPLACEMENT",4); CHKERRQ(ierr);
  ierr = mField.set_field_order(0,MBEDGE,"DISPLACEMENT",4); CHKERRQ(ierr);
  ierr = mField.set_field_order(0,MBVERTEX,"DISPLACEMENT",1); CHKERRQ(ierr);

  //build field
  ierr = mField.build_fields(); CHKERRQ(ierr);
  /*ierr = mField.set_field_order(0,MBTET,"DISPLACEMENT",1); CHKERRQ(ierr);
  ierr = mField.set_field_order(0,MBTRI,"DISPLACEMENT",1); CHKERRQ(ierr);
  ierr = mField.set_field_order(0,MBEDGE,"DISPLACEMENT",1); CHKERRQ(ierr);
  ierr = mField.build_fields(); CHKERRQ(ierr);*/

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

  struct ElasticFEMethod: public FEMethod_UpLevelStudent {

    Mat &Aij;
    Vec& F,Diagonal;
    ElasticFEMethod(
      Interface& _moab,Mat &_Aij,Vec& _F,
      double _lambda,double _mu,Range &_SideSet1,Range &_SideSet2): 
      FEMethod_UpLevelStudent(_moab,1),Aij(_Aij),F(_F),
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

      //VEC & MAT Options
      //If index is set to -1 ingonre its assembly
      VecSetOption(F, VEC_IGNORE_NEGATIVE_INDICES,PETSC_TRUE); 
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

    vector<double> g_NTET,g_NTRI;
    
    PetscErrorCode preProcess() {
      PetscFunctionBegin;
      FEMethod_LowLevelStudent::preProcess();
      PetscSynchronizedPrintf(PETSC_COMM_WORLD,"Start Assembly\n",pcomm->rank(),v2-v1,t2-t1);
      ierr = PetscGetTime(&v1); CHKERRQ(ierr);
      ierr = PetscGetCPUTime(&t1); CHKERRQ(ierr);
      g_NTET.resize(4*45);
      ShapeMBTET(&g_NTET[0],G_TET_X45,G_TET_Y45,G_TET_Z45,45);
      g_NTRI.resize(3*7);
      ShapeMBTRI(&g_NTRI[0],G_TRI_X7,G_TRI_Y7,7); 
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
      for(int rr = 0;rr<6;rr++) {
	D_mu(rr,rr) = rr<3 ? 2 : 1;
      }
      D = lambda*D_lambda + mu*D_mu;
      //cerr << D_lambda << endl;
      //cerr << D_mu << endl;
      //cerr << D << endl;
      ierr = VecDuplicate(F,&Diagonal); CHKERRQ(ierr);
      PetscFunctionReturn(0);
    }

    PetscErrorCode postProcess() {
      PetscFunctionBegin;
      ierr = VecGhostUpdateBegin(F,ADD_VALUES,SCATTER_REVERSE); CHKERRQ(ierr);
      ierr = VecGhostUpdateEnd(F,ADD_VALUES,SCATTER_REVERSE); CHKERRQ(ierr);
      ierr = VecAssemblyBegin(F); CHKERRQ(ierr);
      ierr = VecAssemblyEnd(F); CHKERRQ(ierr);
      ierr = VecAssemblyBegin(Diagonal); CHKERRQ(ierr);
      ierr = VecAssemblyEnd(Diagonal); CHKERRQ(ierr);
      ierr = MatDiagonalSet(Aij,Diagonal,ADD_VALUES); CHKERRQ(ierr);
      ierr = MatAssemblyBegin(Aij,MAT_FINAL_ASSEMBLY); CHKERRQ(ierr);
      ierr = MatAssemblyEnd(Aij,MAT_FINAL_ASSEMBLY); CHKERRQ(ierr);
      ierr = PetscGetTime(&v2); CHKERRQ(ierr);
      ierr = PetscGetCPUTime(&t2); CHKERRQ(ierr);
      ierr = VecDestroy(&Diagonal); CHKERRQ(ierr);
      PetscSynchronizedPrintf(PETSC_COMM_WORLD,"End Assembly: Rank %d Time = %f CPU Time = %f\n",pcomm->rank(),v2-v1,t2-t1);
      PetscFunctionReturn(0);
    }

    PetscErrorCode NeumannBC() {
      PetscFunctionBegin;
      
      ublas::vector<FieldData,ublas::bounded_array<double,3> > traction(3);
      traction[0] = 0;
      traction[1] = 1;
      traction[2] = 0;

      SideNumber_multiIndex& side_table = const_cast<SideNumber_multiIndex&>(fe_ent_ptr->get_side_number_table());
      SideNumber_multiIndex::nth_index<1>::type::iterator siit = side_table.get<1>().lower_bound(boost::make_tuple(MBTRI,0));
      SideNumber_multiIndex::nth_index<1>::type::iterator hi_siit = side_table.get<1>().upper_bound(boost::make_tuple(MBTRI,4));
      for(;siit!=hi_siit;siit++) {
	Range::iterator fit = find(SideSet2.begin(),SideSet2.end(),siit->ent);
	if(fit==SideSet2.end()) continue;

	ierr = ShapeFunctions_TRI(siit->ent,g_NTRI);  CHKERRQ(ierr);

	//nodes
	vector<DofIdx>& RowGlob_nodes = RowGlob[0];
	vector< ublas::matrix<FieldData> > FaceNMatrix_nodes;
	ierr = GetGaussRowFaceNMatrix(siit->ent,"DISPLACEMENT",FaceNMatrix_nodes,MBVERTEX); CHKERRQ(ierr);
	//copy(FaceNMatrix_nodes.begin(),FaceNMatrix_nodes.end(),ostream_iterator<ublas::matrix<FieldData> >(cerr," \n")); cerr << endl;

	//edges
	vector<vector<ublas::matrix<FieldData> > > FaceNMatrix_edges(6);
	SideNumber_multiIndex::nth_index<1>::type::iterator siiit = side_table.get<1>().lower_bound(boost::make_tuple(MBEDGE,0));
	SideNumber_multiIndex::nth_index<1>::type::iterator hi_siiit = side_table.get<1>().upper_bound(boost::make_tuple(MBEDGE,6));
	for(;siiit!=hi_siiit;siiit++) {
	  if(RowGlob[1+siiit->side_number].size()>0) {
	    ierr = GetGaussRowFaceNMatrix(siit->ent,"DISPLACEMENT",FaceNMatrix_edges[siiit->side_number],MBEDGE,siiit->ent); CHKERRQ(ierr);
	    //cerr << "ee ";
	    //copy(FaceNMatrix_face.begin(),FaceNMatrix_face.end(),ostream_iterator<ublas::matrix<FieldData> >(cerr," \n")); cerr << endl;
	  }
	}

	//faces
	vector<DofIdx>& RowGlob_face = RowGlob[1+6+siit->side_number];
	vector< ublas::matrix<FieldData> > FaceNMatrix_face;
	if(RowGlob_face.size()>0) {
	  ierr = GetGaussRowFaceNMatrix(siit->ent,"DISPLACEMENT",FaceNMatrix_face,MBTRI); CHKERRQ(ierr);
	  //copy(FaceNMatrix_face.begin(),FaceNMatrix_face.end(),ostream_iterator<ublas::matrix<FieldData> >(cerr," \n")); cerr << endl;
	}

	//normal and area of trianagular face
	const EntityHandle *conn_face;
	int num_nodes_face;
	rval = moab.get_connectivity(siit->ent,conn_face,num_nodes_face,true); CHKERR_PETSC(rval);
	double coords_face[9];
	rval = moab.get_coords(conn_face,num_nodes_face,coords_face); CHKERR_PETSC(rval);
	ierr = ShapeDiffMBTRI(diffNTRI); CHKERRQ(ierr);
	double normal[3];
	ierr = ShapeFaceNormalMBTRI(diffNTRI,coords_face,normal); CHKERRQ(ierr);
	double area = cblas_dnrm2(3,normal,1);

	//calulate & assemble

	//nodes
	ublas::vector<FieldData> f_ext_nodes = ublas::zero_vector<FieldData>(FaceNMatrix_nodes[0].size2());
	int g_dim = get_dim_gNTRI();
	for(int gg = 0;gg<g_dim;gg++) {
	  double w = area*G_TRI_W7[gg];
	  f_ext_nodes += w*prod(trans(FaceNMatrix_nodes[gg]), traction);
	}
	ierr = VecSetValues(F,RowGlob_nodes.size(),&(RowGlob_nodes)[0],&(f_ext_nodes.data())[0],ADD_VALUES); CHKERRQ(ierr);

	//edges
	siiit = side_table.get<1>().lower_bound(boost::make_tuple(MBEDGE,0));
	for(;siiit!=hi_siiit;siiit++) {
	  vector<DofIdx>& RowGlob_edge = RowGlob[1+siiit->side_number];
	  if(RowGlob_edge.size()>0) {
	    vector<ublas::matrix<FieldData> >& FaceNMatrix_edge = FaceNMatrix_edges[siiit->side_number];
	    ublas::vector<FieldData> f_ext_edges = ublas::zero_vector<FieldData>(FaceNMatrix_edge[0].size2());
	    for(int gg = 0;gg<g_dim;gg++) {
	      double w = area*G_TRI_W7[gg];
	      f_ext_edges += w*prod(trans(FaceNMatrix_edge[gg]), traction);
	    }
	    if(RowGlob_edge.size()!=f_ext_edges.size()) SETERRQ(PETSC_COMM_SELF,1,"wrong size: RowGlob_edge.size()!=f_ext_edges.size()");
	    ierr = VecSetValues(F,RowGlob_edge.size(),&(RowGlob_edge[0]),&(f_ext_edges.data())[0],ADD_VALUES); CHKERRQ(ierr);
	  }
	}
  
 	//face     
	if(RowGlob_face.size()>0) {
	  ublas::vector<FieldData> f_ext_faces = ublas::zero_vector<FieldData>(FaceNMatrix_face[0].size2());
	  for(int gg = 0;gg<g_dim;gg++) {
	    double w = area*G_TRI_W7[gg];
	    f_ext_faces += w*prod(trans(FaceNMatrix_face[gg]), traction);
	  }
	  if(RowGlob_face.size()!=f_ext_faces.size()) SETERRQ(PETSC_COMM_SELF,1,"wrong size: RowGlob_face.size()!=f_ext_faces.size()");
	  ierr = VecSetValues(F,RowGlob_face.size(),&(RowGlob_face)[0],&(f_ext_faces.data())[0],ADD_VALUES); CHKERRQ(ierr);
	}


      }

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
	  //cerr << rowDiffNMatrices[row_mat][0] << endl;
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
	ierr = GetGaussRowNMatrix("DISPLACEMENT",MBTET,rowNMatrices[row_mat]); CHKERRQ(ierr);
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
	ierr = GetGaussColNMatrix("DISPLACEMENT",MBTET,colNMatrices[col_mat]); CHKERRQ(ierr);
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
	    vector<DofIdx>::iterator it = find(ColGlob[cc].begin(),ColGlob[cc].end(),riit->get_petsc_gloabl_dof_idx());
	    if( it!=ColGlob[cc].end() ) *it = -1; // of idx is set -1 column is not assembled
	  }
	  for(int rr = 0;rr<row_mat;rr++) {
	    vector<DofIdx>::iterator it = find(RowGlob[rr].begin(),RowGlob[rr].end(),riit->get_petsc_gloabl_dof_idx());
	    if( it!=RowGlob[rr].end() ) *it = -1; // of idx is set -1 row is not assembled
	  }
	}
      }
      if(DirihletBC.size()>0) {
	DirihletBCDiagVal.resize(DirihletBC.size());
	fill(DirihletBCDiagVal.begin(),DirihletBCDiagVal.end(),1);
	ierr = VecSetValues(Diagonal,DirihletBC.size(),&(DirihletBC[0]),&DirihletBCDiagVal[0],INSERT_VALUES); CHKERRQ(ierr);
      }

      PetscFunctionReturn(0);
    }

    PetscErrorCode operator()() {
      PetscFunctionBegin;
      ierr = OpStudentStart_TET(g_NTET); CHKERRQ(ierr);
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
	  ublas::matrix<FieldData> &row_Mat = (rowNMatrices[rr])[gg];
	  if(gg == 0) f[rr] = ublas::zero_vector<FieldData>(row_Mat.size2());
	  ///f matrices
	  // calulate body force (f.e. garvity force
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

      //Neumann Boundary Conditions
      ierr = NeumannBC(); CHKERRQ(ierr);

      ierr = OpStudentEnd(); CHKERRQ(ierr);
      PetscFunctionReturn(0);
    }

  };

  Range SideSet1,SideSet2;
  ierr = mField.get_Cubit_msId_entities_by_dimension(1,SideSet,2,SideSet1,true); CHKERRQ(ierr);
  ierr = mField.get_Cubit_msId_entities_by_dimension(2,SideSet,2,SideSet2,true); CHKERRQ(ierr);
  PetscPrintf(PETSC_COMM_WORLD,"Nb. faces in SideSet 1 : %u\n",SideSet1.size());
  PetscPrintf(PETSC_COMM_WORLD,"Nb. faces in SideSet 1 : %u\n",SideSet2.size());

  const double YoungModulus = 1;
  const double PoissonRatio = 0.25;
  ElasticFEMethod MyFE(moab,Aij,F,LAMBDA(YoungModulus,PoissonRatio),MU(YoungModulus,PoissonRatio),SideSet1,SideSet2);
  ierr = mField.loop_finite_elements("ELASTIC_MECHANICS","ELASTIC",MyFE);  CHKERRQ(ierr);
  PetscSynchronizedFlush(PETSC_COMM_WORLD);

  //Matrix View
  //MatView(Aij,PETSC_VIEWER_DRAW_WORLD);//PETSC_VIEWER_STDOUT_WORLD);
  //std::string wait;
  //std::cin >> wait;

  //Solver
  KSP solver;
  ierr = KSPCreate(PETSC_COMM_WORLD,&solver); CHKERRQ(ierr);
  ierr = KSPSetOperators(solver,Aij,Aij,SAME_NONZERO_PATTERN); CHKERRQ(ierr);
  ierr = KSPSetFromOptions(solver); CHKERRQ(ierr);
  ierr = KSPSetUp(solver); CHKERRQ(ierr);

  Vec D;
  ierr = VecDuplicate(F,&D); CHKERRQ(ierr);
  ierr = KSPSolve(solver,F,D); CHKERRQ(ierr);
  ierr = VecGhostUpdateBegin(D,INSERT_VALUES,SCATTER_FORWARD); CHKERRQ(ierr);
  ierr = VecGhostUpdateEnd(D,INSERT_VALUES,SCATTER_FORWARD); CHKERRQ(ierr);

  //Save data on mesh
  ierr = mField.set_global_VecCreateGhost("ELASTIC_MECHANICS",Row,D,INSERT_VALUES,SCATTER_REVERSE); CHKERRQ(ierr);
  //ierr = VecView(F,PETSC_VIEWER_STDOUT_WORLD); CHKERRQ(ierr);

  // Write Displacements DOFS on Vertices
  struct MyEntMethod: public moabField::EntMethod {
    ErrorCode rval;
    PetscErrorCode ierr;

    Tag th_disp;
    MyEntMethod(Interface& _moab): EntMethod(_moab) {
      double def_VAL[3] = {0,0,0};
      // create TAG
      rval = moab.tag_get_handle("DISPLACEMENTS_VAL",3,MB_TYPE_DOUBLE,th_disp,MB_TAG_CREAT|MB_TAG_SPARSE,&def_VAL); CHKERR(rval);
    }

    vector<double> vals;
    Range nodes;
    
    PetscErrorCode preProcess() {
      PetscFunctionBegin;
      PetscPrintf(PETSC_COMM_WORLD,"Start postporcess\n");
      rval = moab.get_entities_by_type(0,MBVERTEX,nodes); CHKERR_PETSC(rval);
      vals.resize(nodes.size()*3);
      PetscFunctionReturn(0);
    }
    
    PetscErrorCode operator()() {
      PetscFunctionBegin;
      if(dof_ptr->get_ent_type()!=MBVERTEX) PetscFunctionReturn(0);
      EntityHandle ent = dof_ptr->get_ent();
      int dof_rank = dof_ptr->get_dof_rank();
      double fval = dof_ptr->get_FieldData();
      Range::iterator nit = find(nodes.begin(),nodes.end(),ent);
      if(nit==nodes.end()) SETERRQ(PETSC_COMM_SELF,1,"data inconsistency");
      unsigned int pos = std::distance(nodes.begin(),nit);
      pos = 3*pos+dof_rank;
      if(pos>vals.size()) SETERRQ(PETSC_COMM_SELF,1,"data inconsistency");
      vals[pos] = fval;
      //cerr << pos << " --> " << fval << " ent " << ent << endl;
      PetscFunctionReturn(0);
    }

    PetscErrorCode postProcess() {
      PetscFunctionBegin;
      ierr = moab.tag_set_data(th_disp,nodes,&vals[0]); CHKERRQ(ierr);
      PetscPrintf(PETSC_COMM_WORLD,"End postporcess\n");
      PetscFunctionReturn(0);
    }

  };

  MyEntMethod ent_method(moab);
  ierr = mField.loop_dofs("ELASTIC_MECHANICS","DISPLACEMENT",Row,ent_method); CHKERRQ(ierr);

  if(pcomm->rank()==0) {
    EntityHandle out_meshset;
    rval = moab.create_meshset(MESHSET_SET,out_meshset); CHKERR_PETSC(rval);
    ierr = mField.problem_get_FE("ELASTIC_MECHANICS","ELASTIC",out_meshset); CHKERRQ(ierr);
    rval = moab.write_file("out.vtk","VTK","",&out_meshset,1); CHKERR_PETSC(rval);
    rval = moab.delete_entities(&out_meshset,1); CHKERR_PETSC(rval);
  }

  struct PostProcFEMethod: public FEMethod_UpLevelStudent {

    Interface& moab_post_proc;
    Core mb_instance_post_proc;
    Interface& moab_ref;
    Core mb_instance_ref;

    const int max_level;
    vector<EntityHandle> meshset_level;
    bool init_ref;

    ParallelComm* pcomm;
    PetscLogDouble t1,t2;
    PetscLogDouble v1,v2;

    PostProcFEMethod(Interface& _moab): FEMethod_UpLevelStudent(_moab,1),
      moab_post_proc(mb_instance_post_proc),moab_ref(mb_instance_ref),
      max_level(2),init_ref(false) {
      pcomm = ParallelComm::get_pcomm(&moab,MYPCOMM_INDEX);

      meshset_level.resize(max_level+1);

    }

    Tag th_disp,th_strain;
    vector<double> g_NTET;

    PetscErrorCode preProcess() {
      PetscFunctionBegin;
      PetscSynchronizedPrintf(PETSC_COMM_WORLD,"Start PostProc\n",pcomm->rank(),v2-v1,t2-t1);
      ierr = PetscGetTime(&v1); CHKERRQ(ierr);
      ierr = PetscGetCPUTime(&t1); CHKERRQ(ierr);
      FEMethod_LowLevelStudent::preProcess();

      if(init_ref) PetscFunctionReturn(0);
      
      double base_coords[] = {
	0,0,0,
	1,0,0,
	0,1,0,
	0,0,1 };
      
      EntityHandle nodes[4];
      for(int nn = 0;nn<4;nn++) {
	rval = moab_ref.create_vertex(&base_coords[3*nn],nodes[nn]); CHKERR_PETSC(rval);
      }
      EntityHandle tet;
      rval = moab_ref.create_element(MBTET,nodes,4,tet); CHKERR_PETSC(rval);

      moabField_Core core_ref(moab_ref);
      moabField& mField_ref = core_ref;
      ierr = mField_ref.seed_ref_level_3D(0,BitRefLevel().set(0)); CHKERRQ(ierr);

      for(int ll = 0;ll<max_level;ll++) {
	PetscPrintf(PETSC_COMM_WORLD,"Refine Level %d\n",ll);
	rval = moab_ref.create_meshset(MESHSET_SET,meshset_level[ll]); CHKERR_PETSC(rval);
	ierr = mField_ref.refine_get_ents(BitRefLevel().set(ll),meshset_level[ll]); CHKERRQ(ierr);
	ierr = mField_ref.add_verices_in_the_middel_of_edges(meshset_level[ll],BitRefLevel().set(ll+1)); CHKERRQ(ierr);
	ierr = mField_ref.refine_TET(meshset_level[ll],BitRefLevel().set(ll+1)); CHKERRQ(ierr);
      }
      rval = moab_ref.create_meshset(MESHSET_SET,meshset_level[max_level]); CHKERR_PETSC(rval);
      ierr = mField_ref.refine_get_ents(BitRefLevel().set(max_level),meshset_level[max_level]); CHKERRQ(ierr);

      std::vector<double> ref_coords;
      rval = moab_ref.get_vertex_coordinates(ref_coords); CHKERR_PETSC(rval);
      g_NTET.resize(4*ref_coords.size()/3);
      ShapeMBTET(&g_NTET[0],&ref_coords[0],&ref_coords[ref_coords.size()/3],&ref_coords[2*ref_coords.size()/3],ref_coords.size()/3);

      double def_VAL[9] = {0,0,0, 0,0,0, 0,0,0};
      // create TAG
      rval = moab_post_proc.tag_get_handle("DISPLACEMENTS_VAL",3,MB_TYPE_DOUBLE,th_disp,MB_TAG_CREAT|MB_TAG_SPARSE,&def_VAL); CHKERR(rval);
      rval = moab_post_proc.tag_get_handle("STRAIN_VAL",9,MB_TYPE_DOUBLE,th_strain,MB_TAG_CREAT|MB_TAG_SPARSE,&def_VAL); CHKERR(rval);

      init_ref = true;

      PetscFunctionReturn(0);
    }

    PetscErrorCode operator()() {
      PetscFunctionBegin;
      ierr = OpStudentStart_TET(g_NTET); CHKERRQ(ierr);

      Range ref_nodes;
      rval = moab_ref.get_entities_by_type(meshset_level[max_level],MBVERTEX,ref_nodes); CHKERR_PETSC(rval);
      if(4*ref_nodes.size()!=g_NTET.size()) SETERRQ(PETSC_COMM_SELF,1,"data inconsistency");
      if(ref_nodes.size()!=coords_at_Gauss_nodes.size()) SETERRQ(PETSC_COMM_SELF,1,"data inconsistency");
      map<EntityHandle,EntityHandle> node_map;
      Range::iterator nit = ref_nodes.begin();
      for(int nn = 0;nit!=ref_nodes.end();nit++,nn++) {
	EntityHandle &node = node_map[*nit];
	rval = moab_post_proc.create_vertex(&(coords_at_Gauss_nodes[nn]).data()[0],node); CHKERR_PETSC(rval);
      }
      Range ref_tets;
      rval = moab_ref.get_entities_by_type(meshset_level[max_level],MBTET,ref_tets); CHKERR_PETSC(rval);
      Range::iterator tit = ref_tets.begin();
      for(;tit!=ref_tets.end();tit++) {
	const EntityHandle *conn_ref;
        int num_nodes;
	rval = moab_ref.get_connectivity(*tit,conn_ref,num_nodes,true); CHKERR_PETSC(rval);
	EntityHandle conn_post_proc[num_nodes];
	for(int nn = 0;nn<num_nodes;nn++) {
	  conn_post_proc[nn] = node_map[conn_ref[nn]];
	}
	EntityHandle ref_tet;
	rval = moab_post_proc.create_element(MBTET,conn_post_proc,4,ref_tet); CHKERR_PETSC(rval);
      }

      Data_at_Gauss_pt::iterator diit = data_at_gauss_pt.find("DISPLACEMENT");
      if(diit==data_at_gauss_pt.end()) SETERRQ(PETSC_COMM_SELF,1,"no DISPLACEMENT !!!");
      vector< ublas::vector<FieldData> > &data = diit->second;
      vector< ublas::vector<FieldData> >::iterator vit = data.begin();
      map<EntityHandle,EntityHandle>::iterator mit = node_map.begin();
      for(;vit!=data.end();vit++,mit++) {
	rval = moab_post_proc.tag_set_data(th_disp,&mit->second,1,&vit->data()[0]); CHKERR_PETSC(rval);
      }

      //Strains to Noades in PostProc Mesh
      vector< ublas::matrix< FieldData > > GradU_at_GaussPt;
      ierr = GetGaussDiffDataVector("DISPLACEMENT",GradU_at_GaussPt); CHKERRQ(ierr);
      vector< ublas::matrix< FieldData > >::iterator viit = GradU_at_GaussPt.begin();
      mit = node_map.begin();
      for(;viit!=GradU_at_GaussPt.end();viit++,mit++) {
	ublas::matrix< FieldData > Strain = 0.5*( (*viit) + trans(*viit) );
	rval = moab_post_proc.tag_set_data(th_strain,&mit->second,1,&(Strain.data()[0])); CHKERR_PETSC(rval);
      }

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

  PostProcFEMethod fe_post_proc_method(moab);
  ierr = mField.loop_finite_elements("ELASTIC_MECHANICS","ELASTIC",fe_post_proc_method);  CHKERRQ(ierr);
  PetscSynchronizedFlush(PETSC_COMM_WORLD);
  if(pcomm->rank()==0) {
    rval = fe_post_proc_method.moab_post_proc.write_file("out_post_proc.vtk","VTK",""); CHKERR_PETSC(rval);
  }

  //detroy matrices
  ierr = VecDestroy(&F); CHKERRQ(ierr);
  ierr = VecDestroy(&D); CHKERRQ(ierr);
  ierr = MatDestroy(&Aij); CHKERRQ(ierr);
  ierr = KSPDestroy(&solver); CHKERRQ(ierr);


  ierr = PetscGetTime(&v2);CHKERRQ(ierr);
  ierr = PetscGetCPUTime(&t2);CHKERRQ(ierr);

  PetscSynchronizedPrintf(PETSC_COMM_WORLD,"Total Rank %d Time = %f CPU Time = %f\n",pcomm->rank(),v2-v1,t2-t1);
  PetscSynchronizedFlush(PETSC_COMM_WORLD);

  PetscFinalize();

}

