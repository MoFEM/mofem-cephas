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

#ifndef __ELASTICFEMETHOD_HPP__
#define __ELASTICFEMETHOD_HPP__

namespace MoFEM {

struct ElasticFEMethod: public FEMethod_UpLevelStudent {

    Mat Aij;
    Vec F,Diagonal;
    Range dummy;

    ElasticFEMethod(
      Interface& _moab): FEMethod_UpLevelStudent(_moab,1),
      Aij(PETSC_NULL),F(PETSC_NULL),SideSet1(dummy),SideSet2(dummy) {};

    ElasticFEMethod(
      Interface& _moab,Mat &_Aij,Vec& _F,
      double _lambda,double _mu,Range &_SideSet1,Range &_SideSet2): 
      FEMethod_UpLevelStudent(_moab,1),Aij(_Aij),F(_F),Diagonal(PETSC_NULL),
      lambda(_lambda),mu(_mu),
      SideSet1(_SideSet1),SideSet2(_SideSet2) { 
      pcomm = ParallelComm::get_pcomm(&moab,MYPCOMM_INDEX);

      RowGlob.resize(1+6+4+1);
      RowLocal.resize(1+6+4+1);
      rowNMatrices.resize(1+6+4+1);
      rowDiffNMatrices.resize(1+6+4+1);
      rowBMatrices.resize(1+6+4+1);
      ColGlob.resize(1+6+4+1);
      ColLocal.resize(1+6+4+1);
      colNMatrices.resize(1+6+4+1);
      colDiffNMatrices.resize(1+6+4+1);
      colBMatrices.resize(1+6+4+1);

      if(F!=PETSC_NULL) {
	//
	Range SideSet1Edges,SideSet1Nodes;
	rval = moab.get_adjacencies(SideSet1,1,false,SideSet1Edges,Interface::UNION); CHKERR_THROW(rval);
	rval = moab.get_connectivity(SideSet1,SideSet1Nodes,true); CHKERR_THROW(rval);
	SideSet1_.insert(SideSet1.begin(),SideSet1.end());
	SideSet1_.insert(SideSet1Edges.begin(),SideSet1Edges.end());
	SideSet1_.insert(SideSet1Nodes.begin(),SideSet1Nodes.end());

	//VEC & MAT Options
	//If index is set to -1 ingonre its assembly
	VecSetOption(F, VEC_IGNORE_NEGATIVE_INDICES,PETSC_TRUE); 
      }

    }; 

    ErrorCode rval;
    
    ParallelComm* pcomm;
    PetscLogDouble t1,t2;
    PetscLogDouble v1,v2;

    double lambda,mu;
    ublas::matrix<FieldData> D_lambda,D_mu,D;

    Range& SideSet1;
    Range& SideSet2;
    Range SideSet1_;

    int row_mat,col_mat;
    vector<vector<DofIdx> > RowGlob;
    vector<vector<DofIdx> > RowLocal;
    vector<vector<ublas::matrix<FieldData> > > rowNMatrices;
    vector<vector<ublas::matrix<FieldData> > > rowDiffNMatrices;
    vector<vector<ublas::matrix<FieldData> > > rowBMatrices;
    vector<vector<DofIdx> > ColGlob;
    vector<vector<DofIdx> > ColLocal;
    vector<vector<ublas::matrix<FieldData> > > colNMatrices;
    vector<vector<ublas::matrix<FieldData> > > colDiffNMatrices;
    vector<vector<ublas::matrix<FieldData> > > colBMatrices;

    vector<DofIdx> DirihletBC;
    vector<FieldData> DirihletBCDiagVal;

    vector<double> g_NTET,g_NTRI;
    
    PetscErrorCode preProcess() {
      PetscFunctionBegin;
      PetscSynchronizedPrintf(PETSC_COMM_WORLD,"Start Assembly\n");
      ierr = PetscGetTime(&v1); CHKERRQ(ierr);
      ierr = PetscGetCPUTime(&t1); CHKERRQ(ierr);
      g_NTET.resize(4*45);
      ShapeMBTET(&g_NTET[0],G_TET_X45,G_TET_Y45,G_TET_Z45,45);
      g_NTRI.resize(3*13);
      ShapeMBTRI(&g_NTRI[0],G_TRI_X13,G_TRI_Y13,13); 
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
      if(Diagonal!=PETSC_NULL) {
	ierr = VecAssemblyBegin(Diagonal); CHKERRQ(ierr);
	ierr = VecAssemblyEnd(Diagonal); CHKERRQ(ierr);
	ierr = MatDiagonalSet(Aij,Diagonal,ADD_VALUES); CHKERRQ(ierr);
	ierr = VecDestroy(&Diagonal); CHKERRQ(ierr);
      }
      ierr = MatAssemblyBegin(Aij,MAT_FINAL_ASSEMBLY); CHKERRQ(ierr);
      ierr = MatAssemblyEnd(Aij,MAT_FINAL_ASSEMBLY); CHKERRQ(ierr);
      ierr = PetscGetTime(&v2); CHKERRQ(ierr);
      ierr = PetscGetCPUTime(&t2); CHKERRQ(ierr);
      PetscSynchronizedPrintf(PETSC_COMM_WORLD,"End Assembly: Rank %d Time = %f CPU Time = %f\n",pcomm->rank(),v2-v1,t2-t1);
      PetscFunctionReturn(0);
    }

    PetscErrorCode NeumannBC(Vec F_ext,ublas::vector<FieldData,ublas::bounded_array<double,3> > &traction,Range& SideSet) {
      PetscFunctionBegin;

      SideNumber_multiIndex& side_table = const_cast<SideNumber_multiIndex&>(fe_ent_ptr->get_side_number_table());
      SideNumber_multiIndex::nth_index<1>::type::iterator siit = side_table.get<1>().lower_bound(boost::make_tuple(MBTRI,0));
      SideNumber_multiIndex::nth_index<1>::type::iterator hi_siit = side_table.get<1>().upper_bound(boost::make_tuple(MBTRI,4));
      for(;siit!=hi_siit;siit++) {
	Range::iterator fit = find(SideSet.begin(),SideSet.end(),siit->ent);
	if(fit==SideSet.end()) continue;

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
	  double w = area*G_TRI_W13[gg];
	  f_ext_nodes += w*prod(trans(FaceNMatrix_nodes[gg]), traction);
	}
	ierr = VecSetValues(F_ext,RowGlob_nodes.size(),&(RowGlob_nodes)[0],&(f_ext_nodes.data())[0],ADD_VALUES); CHKERRQ(ierr);

	//edges
	siiit = side_table.get<1>().lower_bound(boost::make_tuple(MBEDGE,0));
	for(;siiit!=hi_siiit;siiit++) {
	  vector<DofIdx>& RowGlob_edge = RowGlob[1+siiit->side_number];
	  if(RowGlob_edge.size()>0) {
	    vector<ublas::matrix<FieldData> >& FaceNMatrix_edge = FaceNMatrix_edges[siiit->side_number];
	    ublas::vector<FieldData> f_ext_edges = ublas::zero_vector<FieldData>(FaceNMatrix_edge[0].size2());
	    for(int gg = 0;gg<g_dim;gg++) {
	      double w = area*G_TRI_W13[gg];
	      f_ext_edges += w*prod(trans(FaceNMatrix_edge[gg]), traction);
	    }
	    if(RowGlob_edge.size()!=f_ext_edges.size()) SETERRQ(PETSC_COMM_SELF,1,"wrong size: RowGlob_edge.size()!=f_ext_edges.size()");
	    ierr = VecSetValues(F_ext,RowGlob_edge.size(),&(RowGlob_edge[0]),&(f_ext_edges.data())[0],ADD_VALUES); CHKERRQ(ierr);
	  }
	}
  
 	//face     
	if(RowGlob_face.size()>0) {
	  ublas::vector<FieldData> f_ext_faces = ublas::zero_vector<FieldData>(FaceNMatrix_face[0].size2());
	  for(int gg = 0;gg<g_dim;gg++) {
	    double w = area*G_TRI_W13[gg];
	    f_ext_faces += w*prod(trans(FaceNMatrix_face[gg]), traction);
	  }
	  if(RowGlob_face.size()!=f_ext_faces.size()) SETERRQ(PETSC_COMM_SELF,1,"wrong size: RowGlob_face.size()!=f_ext_faces.size()");
	  ierr = VecSetValues(F_ext,RowGlob_face.size(),&(RowGlob_face)[0],&(f_ext_faces.data())[0],ADD_VALUES); CHKERRQ(ierr);
	}


      }

      PetscFunctionReturn(0);

    }

    PetscErrorCode NeumannBC() {
      PetscFunctionBegin;
      
      ublas::vector<FieldData,ublas::bounded_array<double,3> > traction(3);
      traction[0] = 0;
      traction[1] = 1;
      traction[2] = 0;

      ierr = NeumannBC(F,traction,SideSet2); CHKERRQ(ierr);

      PetscFunctionReturn(0);
    }

    PetscErrorCode ApplyDirihletBC() {
      PetscFunctionBegin;
      //Dirihlet form SideSet1
      DirihletBC.resize(0);
      Range::iterator siit1 = SideSet1_.begin();
      for(;siit1!=SideSet1_.end();siit1++) {
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
      PetscFunctionReturn(0);
    }
    
    PetscErrorCode GetMatrices() {
      PetscFunctionBegin;
      //indicies ROWS
      row_mat = 0;
      ierr = GetRowGlobalIndices("DISPLACEMENT",RowGlob[row_mat]); CHKERRQ(ierr);
      ierr = GetRowLocalIndices("DISPLACEMENT",RowLocal[row_mat]); CHKERRQ(ierr);
      ierr = GetGaussRowNMatrix("DISPLACEMENT",rowNMatrices[row_mat]); CHKERRQ(ierr);
      ierr = GetGaussRowDiffNMatrix("DISPLACEMENT",rowDiffNMatrices[row_mat]); CHKERRQ(ierr);
      ierr = MakeBMatrix3D("DISPLACEMENT",rowDiffNMatrices[row_mat],rowBMatrices[row_mat]);  CHKERRQ(ierr);
      row_mat++;
      for(int ee = 0;ee<6;ee++) { //edges matrices
	ierr = GetRowGlobalIndices("DISPLACEMENT",MBEDGE,RowGlob[row_mat],ee); CHKERRQ(ierr);
	ierr = GetRowLocalIndices("DISPLACEMENT",MBEDGE,RowLocal[row_mat],ee); CHKERRQ(ierr);
	if(RowGlob[row_mat].size()!=0) {
	  ierr = GetGaussRowNMatrix("DISPLACEMENT",MBEDGE,rowNMatrices[row_mat],ee); CHKERRQ(ierr);
	  ierr = GetGaussRowDiffNMatrix("DISPLACEMENT",MBEDGE,rowDiffNMatrices[row_mat],ee); CHKERRQ(ierr);
	  ierr = MakeBMatrix3D("DISPLACEMENT",rowDiffNMatrices[row_mat],rowBMatrices[row_mat]);  CHKERRQ(ierr);
	  //cerr << rowDiffNMatrices[row_mat][0] << endl;
	  row_mat++;
	}
      }
      for(int ff = 0;ff<4;ff++) { //faces matrices
	ierr = GetRowGlobalIndices("DISPLACEMENT",MBTRI,RowGlob[row_mat],ff); CHKERRQ(ierr);
	ierr = GetRowLocalIndices("DISPLACEMENT",MBTRI,RowLocal[row_mat],ff); CHKERRQ(ierr);
	if(RowGlob[row_mat].size()!=0) {
	  ierr = GetGaussRowNMatrix("DISPLACEMENT",MBTRI,rowNMatrices[row_mat],ff); CHKERRQ(ierr);
	  ierr = GetGaussRowDiffNMatrix("DISPLACEMENT",MBTRI,rowDiffNMatrices[row_mat],ff); CHKERRQ(ierr);
	  ierr = MakeBMatrix3D("DISPLACEMENT",rowDiffNMatrices[row_mat],rowBMatrices[row_mat]);  CHKERRQ(ierr);
	  row_mat++;
	}
      }
      ierr = GetRowGlobalIndices("DISPLACEMENT",MBTET,RowGlob[row_mat]); CHKERRQ(ierr);
      ierr = GetRowLocalIndices("DISPLACEMENT",MBTET,RowLocal[row_mat]); CHKERRQ(ierr);
      if(RowGlob[row_mat].size() != 0) { //volume matrices
	ierr = GetGaussRowNMatrix("DISPLACEMENT",MBTET,rowNMatrices[row_mat]); CHKERRQ(ierr);
	ierr = GetGaussRowDiffNMatrix("DISPLACEMENT",MBTET,rowDiffNMatrices[row_mat]); CHKERRQ(ierr);
	ierr = MakeBMatrix3D("DISPLACEMENT",rowDiffNMatrices[row_mat],rowBMatrices[row_mat]);  CHKERRQ(ierr);
	row_mat++;
      }

      //indicies COLS
      col_mat = 0;
      ierr = GetColGlobalIndices("DISPLACEMENT",ColGlob[col_mat]); CHKERRQ(ierr);
      ierr = GetColLocalIndices("DISPLACEMENT",ColLocal[col_mat]); CHKERRQ(ierr);
      ierr = GetGaussColNMatrix("DISPLACEMENT",colNMatrices[col_mat]); CHKERRQ(ierr);
      ierr = GetGaussColDiffNMatrix("DISPLACEMENT",colDiffNMatrices[col_mat]); CHKERRQ(ierr);
      ierr = MakeBMatrix3D("DISPLACEMENT",colDiffNMatrices[col_mat],colBMatrices[col_mat]);  CHKERRQ(ierr);
      col_mat++;
      for(int ee = 0;ee<6;ee++) { //edges matrices
	ierr = GetColGlobalIndices("DISPLACEMENT",MBEDGE,ColGlob[col_mat],ee); CHKERRQ(ierr);
	ierr = GetColLocalIndices("DISPLACEMENT",MBEDGE,ColLocal[col_mat],ee); CHKERRQ(ierr);
	if(ColGlob[col_mat].size()!=0) {
	  ierr = GetGaussColNMatrix("DISPLACEMENT",MBEDGE,colNMatrices[col_mat],ee); CHKERRQ(ierr);
	  ierr = GetGaussColDiffNMatrix("DISPLACEMENT",MBEDGE,colDiffNMatrices[col_mat],ee); CHKERRQ(ierr);
	  ierr = MakeBMatrix3D("DISPLACEMENT",colDiffNMatrices[col_mat],colBMatrices[col_mat]);  CHKERRQ(ierr);
	  col_mat++;
	}
      }
      for(int ff = 0;ff<4;ff++) { //faces matrices
	ierr = GetColGlobalIndices("DISPLACEMENT",MBTRI,ColGlob[col_mat],ff); CHKERRQ(ierr);
	ierr = GetColLocalIndices("DISPLACEMENT",MBTRI,ColLocal[col_mat],ff); CHKERRQ(ierr);
	if(ColGlob[col_mat].size()!=0) {
	  ierr = GetGaussColNMatrix("DISPLACEMENT",MBTRI,colNMatrices[col_mat],ff); CHKERRQ(ierr);
	  ierr = GetGaussColDiffNMatrix("DISPLACEMENT",MBTRI,colDiffNMatrices[col_mat],ff); CHKERRQ(ierr);
	  ierr = MakeBMatrix3D("DISPLACEMENT",colDiffNMatrices[col_mat],colBMatrices[col_mat]);  CHKERRQ(ierr);
	  col_mat++;
	}
      }
      ierr = GetColGlobalIndices("DISPLACEMENT",MBTET,ColGlob[col_mat]); CHKERRQ(ierr);
      ierr = GetColLocalIndices("DISPLACEMENT",MBTET,ColLocal[col_mat]); CHKERRQ(ierr);
      if(ColGlob[col_mat].size() != 0) { //volume matrices
	ierr = GetGaussColNMatrix("DISPLACEMENT",MBTET,colNMatrices[col_mat]); CHKERRQ(ierr);
	ierr = GetGaussColDiffNMatrix("DISPLACEMENT",MBTET,colDiffNMatrices[col_mat]); CHKERRQ(ierr);
	ierr = MakeBMatrix3D("DISPLACEMENT",colDiffNMatrices[col_mat],colBMatrices[col_mat]);  CHKERRQ(ierr);
	col_mat++;
      }

      PetscFunctionReturn(0);
    }

    ublas::matrix<ublas::matrix<FieldData> > K;
    PetscErrorCode Stiffness() {
      PetscFunctionBegin;
      K.resize(row_mat,col_mat);
      int g_dim = g_NTET.size()/4;
      for(int rr = 0;rr<row_mat;rr++) {
	for(int cc = 0;cc<col_mat;cc++) {
	  for(int gg = 0;gg<g_dim;gg++) {
	    ublas::matrix<FieldData> &row_Mat = (rowBMatrices[rr])[gg];
	    ublas::matrix<FieldData> &col_Mat = (colBMatrices[cc])[gg];
	    ///K matrices
	    double w = V*G_TET_W45[gg];
	    ublas::matrix<FieldData> BTD = prod( trans(row_Mat), w*D );
	    if(gg == 0) {
	      K(rr,cc) = prod(BTD , col_Mat ); // int BT*D*B
	    } else {
	      K(rr,cc) += prod(BTD , col_Mat ); // int BT*D*B
	    }
	  }
	}
      }
      PetscFunctionReturn(0);
    }

    PetscErrorCode Lhs() {
      PetscFunctionBegin;
      ierr = Stiffness(); CHKERRQ(ierr);
      for(int rr = 0;rr<row_mat;rr++) {
	for(int cc = 0;cc<col_mat;cc++) {
	  if(ColGlob[cc].size()==0) continue;
	  if(RowGlob[rr].size()!=K(rr,cc).size1()) SETERRQ(PETSC_COMM_SELF,1,"data inconsistency");
	  if(ColGlob[cc].size()!=K(rr,cc).size2()) SETERRQ(PETSC_COMM_SELF,1,"data inconsistency");
	  ierr = MatSetValues(Aij,RowGlob[rr].size(),&(RowGlob[rr])[0],ColGlob[cc].size(),&(ColGlob[cc])[0],&(K(rr,cc).data())[0],ADD_VALUES); CHKERRQ(ierr);
	}
      }
      PetscFunctionReturn(0);
    }

    PetscErrorCode Rhs() {
      PetscFunctionBegin;
      ublas::vector<FieldData> f[row_mat];
      int g_dim = g_NTET.size()/4;
      for(int rr = 0;rr<row_mat;rr++) {
	for(int gg = 0;gg<g_dim;gg++) {
	  ublas::matrix<FieldData> &row_Mat = (rowNMatrices[rr])[gg];
	  if(gg == 0) f[rr] = ublas::zero_vector<FieldData>(row_Mat.size2());
	  ///f matrices
	  // calulate body force (f.e. garvity force
	}
	if(RowGlob[rr].size()!=f[rr].size()) SETERRQ(PETSC_COMM_SELF,1,"data inconsistency");
	if(RowGlob[rr].size()==0) continue;
	ierr = VecSetValues(F,RowGlob[rr].size(),&(RowGlob[rr])[0],&(f[rr].data())[0],ADD_VALUES); CHKERRQ(ierr);
      }
      PetscFunctionReturn(0);
    }

    PetscErrorCode Fint() {
      PetscFunctionBegin;
      //Gradient at Gauss points; 
      vector< ublas::matrix< FieldData > > GradU_at_GaussPt;
      ierr = GetGaussDiffDataVector("DISPLACEMENT",GradU_at_GaussPt); CHKERRQ(ierr);
      unsigned int g_dim = g_NTET.size()/4;
      assert(GradU_at_GaussPt.size() == g_dim);
      NOT_USED(g_dim);
      vector< ublas::matrix< FieldData > >::iterator viit = GradU_at_GaussPt.begin();
      int gg = 0;
      for(;viit!=GradU_at_GaussPt.end();viit++,gg++) {
	  ublas::matrix< FieldData > GradU = *viit;
	  ublas::matrix< FieldData > Strain = 0.5*( GradU + trans(GradU) );
	  ublas::vector< FieldData > VoightStrain(6);
	  VoightStrain[0] = Strain(0,0);
	  VoightStrain[1] = Strain(1,1);
	  VoightStrain[2] = Strain(2,2);
	  VoightStrain[3] = 2*Strain(0,1);
	  VoightStrain[4] = 2*Strain(1,2);
	  VoightStrain[5] = 2*Strain(2,0);
	  double w = V*G_TET_W45[gg];
	  ublas::vector<FieldData> VoightStress = prod(w*D,VoightStrain);
	  //BT * VoigtStress
	  for(int rr = 0;rr<row_mat;rr++) {
	    ublas::matrix<FieldData> &B = (rowBMatrices[rr])[gg];
	    ublas::vector<FieldData> f_int = prod( trans(B), VoightStress );
	    if(RowGlob[rr].size()!=f_int.size()) SETERRQ(PETSC_COMM_SELF,1,"data inconsistency");
	    if(RowGlob[rr].size()==0) continue;
	    ierr = VecSetValues(F,RowGlob[rr].size(),&(RowGlob[rr])[0],&(f_int.data()[0]),ADD_VALUES); CHKERRQ(ierr);
	  }
      }

      PetscFunctionReturn(0);
    }

    PetscErrorCode RhsAndLhs() {
      PetscFunctionBegin;
      ierr = OpStudentStart_TET(g_NTET); CHKERRQ(ierr);
      ierr = GetMatrices(); CHKERRQ(ierr);

      ierr = Rhs(); CHKERRQ(ierr);
      ierr = Lhs(); CHKERRQ(ierr);

      ierr = OpStudentEnd(); CHKERRQ(ierr);
      PetscFunctionReturn(0);
    }

    PetscErrorCode operator()() {
      PetscFunctionBegin;
      ierr = OpStudentStart_TET(g_NTET); CHKERRQ(ierr);
      ierr = GetMatrices(); CHKERRQ(ierr);

      //Dirihlet Boundary Condition
      ApplyDirihletBC();
      if(Diagonal!=PETSC_NULL) {
	if(DirihletBC.size()>0) {
	  DirihletBCDiagVal.resize(DirihletBC.size());
	  fill(DirihletBCDiagVal.begin(),DirihletBCDiagVal.end(),1);
	  ierr = VecSetValues(Diagonal,DirihletBC.size(),&(DirihletBC[0]),&DirihletBCDiagVal[0],INSERT_VALUES); CHKERRQ(ierr);
	}
      }

      //Assembly Aij and F
      ierr = RhsAndLhs(); CHKERRQ(ierr);

      //Neumann Boundary Conditions
      ierr = NeumannBC(); CHKERRQ(ierr);

      ierr = OpStudentEnd(); CHKERRQ(ierr);
      PetscFunctionReturn(0);
    }

  };

}

#endif //__ELASTICFEMETHOD_HPP__
