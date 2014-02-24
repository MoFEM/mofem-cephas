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

#include <boost/numeric/ublas/symmetric.hpp>

namespace MoFEM {

struct ElasticFEMethod: public FEMethod_UpLevelStudent {

    FieldInterface& mField;
    Mat Aij;
    Vec Data,F;

    ElasticFEMethod(
      FieldInterface& _mField): FEMethod_UpLevelStudent(_mField.get_moab(),1), mField(_mField),
      Aij(PETSC_NULL),Data(PETSC_NULL),F(PETSC_NULL) {};

    bool propeties_from_BlockSet_Mat_ElasticSet;
    ElasticFEMethod(
      FieldInterface& _mField,BaseDirihletBC *_dirihlet_ptr,Mat &_Aij,Vec &_D,Vec& _F,
      double _lambda,double _mu): 
      FEMethod_UpLevelStudent(_mField.get_moab(),_dirihlet_ptr,1), mField(_mField),
      Aij(_Aij),Data(_D),F(_F),
      lambda(_lambda),mu(_mu) { 
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
	//VEC & MAT Options
	//If index is set to -1 ingonre its assembly
	VecSetOption(F, VEC_IGNORE_NEGATIVE_INDICES,PETSC_TRUE); 
      }

      g_NTET.resize(4*45);
      ShapeMBTET(&g_NTET[0],G_TET_X45,G_TET_Y45,G_TET_Z45,45);
      G_W_TET = G_TET_W45;
      g_NTRI.resize(3*13);
      ShapeMBTRI(&g_NTRI[0],G_TRI_X13,G_TRI_Y13,13); 
      G_W_TRI = G_TRI_W13;

      propeties_from_BlockSet_Mat_ElasticSet = false;
      for(_IT_CUBITMESHSETS_BY_BCDATA_TYPE_FOR_LOOP_(mField,BlockSet|Mat_ElasticSet,it)) {
	propeties_from_BlockSet_Mat_ElasticSet = true;
      }

    }; 

    ErrorCode rval;
    
    ParallelComm* pcomm;
    PetscLogDouble t1,t2;
    PetscLogDouble v1,v2;

    double lambda,mu;
    ublas::matrix<FieldData> D_lambda,D_mu,D;

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

    vector< ublas::matrix<FieldData> > invH;
    vector< FieldData > detH;

    vector<DofIdx> DirihletBC;

    vector<double> g_NTET,g_NTRI;
    const double* G_W_TET;
    const double* G_W_TRI;

    virtual PetscErrorCode calulateD(double _lambda,double _mu) {
      PetscFunctionBegin;

      D = _lambda*D_lambda + _mu*D_mu;
      //cerr << D_lambda << endl;
      //cerr << D_mu << endl;
      //cerr << D << endl;

      PetscFunctionReturn(0);
    }

    PetscErrorCode GetMatParameters(double *_lambda,double *_mu) {
      PetscFunctionBegin;

      *_lambda = lambda;
      *_mu = mu;


      if(propeties_from_BlockSet_Mat_ElasticSet) {
	EntityHandle ent = fe_ptr->get_ent();
	for(_IT_CUBITMESHSETS_BY_BCDATA_TYPE_FOR_LOOP_(mField,BlockSet|Mat_ElasticSet,it)) {

	  Mat_Elastic mydata;
	  ierr = it->get_attribute_data_structure(mydata); CHKERRQ(ierr);

	  Range meshsets;
	  rval = moab.get_entities_by_type(it->meshset,MBENTITYSET,meshsets,true); CHKERR_PETSC(rval);
	  meshsets.insert(it->meshset);
	  for(Range::iterator mit = meshsets.begin();mit != meshsets.end(); mit++) {
	    if( moab.contains_entities(*mit,&ent,1) ) {
	      *_lambda = LAMBDA(mydata.data.Young,mydata.data.Poisson);
	      *_mu = MU(mydata.data.Young,mydata.data.Poisson);
	    PetscFunctionReturn(0);  
	  }
	}

      }

      SETERRQ(PETSC_COMM_SELF,1,
	"Element is not in elestic block, however you run linear elastic analysis with that element\n"
	"top tip: check if you update block sets after mesh refinments or interface insertion");

      }

      PetscFunctionReturn(0);
    }

    
    PetscErrorCode preProcess() {
      PetscFunctionBegin;
      PetscSynchronizedPrintf(PETSC_COMM_WORLD,"Start Assembly\n");
      PetscSynchronizedFlush(PETSC_COMM_WORLD); 
      ierr = PetscTime(&v1); CHKERRQ(ierr);
      ierr = PetscGetCPUTime(&t1); CHKERRQ(ierr);
      g_NTET.resize(4*45);
      ShapeMBTET(&g_NTET[0],G_TET_X45,G_TET_Y45,G_TET_Z45,45);
      G_W_TET = G_TET_W45;
      g_NTRI.resize(3*13);
      ShapeMBTRI(&g_NTRI[0],G_TRI_X13,G_TRI_Y13,13); 
      G_W_TRI = G_TRI_W13;

      // See FEAP - - A Finite Element Analysis Program
      D_lambda.resize(6,6);
      D_lambda.clear();
      for(int rr = 0;rr<3;rr++) {
	for(int cc = 0;cc<3;cc++) {
	  D_lambda(rr,cc) = 1;
	}
      }
      D_mu.resize(6,6);
      D_mu.clear();
      for(int rr = 0;rr<6;rr++) {
	D_mu(rr,rr) = rr<3 ? 2 : 1;
      }

      ierr = VecZeroEntries(Data); CHKERRQ(ierr);
      ierr = dirihlet_bc_method_ptr->SetDirihletBC_to_FieldData(this,Data); CHKERRQ(ierr);

      PetscFunctionReturn(0);
    }

    PetscErrorCode postProcess() {
      PetscFunctionBegin;
      // Note MAT_FLUSH_ASSEMBLY
      ierr = MatAssemblyBegin(Aij,MAT_FLUSH_ASSEMBLY); CHKERRQ(ierr);
      ierr = MatAssemblyEnd(Aij,MAT_FLUSH_ASSEMBLY); CHKERRQ(ierr);
      ierr = dirihlet_bc_method_ptr->SetDirihletBC_to_MatrixDiagonal(this,Aij); CHKERRQ(ierr);
      ierr = VecAssemblyBegin(F); CHKERRQ(ierr);
      ierr = VecAssemblyEnd(F); CHKERRQ(ierr);
      ierr = dirihlet_bc_method_ptr->SetDirihletBC_to_RHS(this,F); CHKERRQ(ierr);
      ierr = VecAssemblyBegin(F); CHKERRQ(ierr);
      ierr = VecAssemblyEnd(F); CHKERRQ(ierr);
      ierr = PetscTime(&v2); CHKERRQ(ierr);
      ierr = PetscGetCPUTime(&t2); CHKERRQ(ierr);
      PetscSynchronizedPrintf(PETSC_COMM_WORLD,"End Assembly: Rank %d Time = %f CPU Time = %f\n",pcomm->rank(),v2-v1,t2-t1);
      PetscSynchronizedFlush(PETSC_COMM_WORLD); 
      PetscFunctionReturn(0);
    }

    virtual PetscErrorCode NeumannBC_Faces(Vec F_ext,
      double pressure,ublas::vector<FieldData,ublas::bounded_array<double,3> > traction,
      Range& faces) {
      PetscFunctionBegin;

      SideNumber_multiIndex& side_table = const_cast<SideNumber_multiIndex&>(fe_ent_ptr->get_side_number_table());
      SideNumber_multiIndex::nth_index<1>::type::iterator siit = side_table.get<1>().lower_bound(boost::make_tuple(MBTRI,0));
      SideNumber_multiIndex::nth_index<1>::type::iterator hi_siit = side_table.get<1>().upper_bound(boost::make_tuple(MBTRI,4));
      for(;siit!=hi_siit;siit++) {
	Range::iterator fit = find(faces.begin(),faces.end(),siit->ent);
	if(fit==faces.end()) continue;

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
	    //copy(FaceNMatrix_edges[siiit->side_number].begin(),FaceNMatrix_edges[siiit->side_number].end(),ostream_iterator<ublas::matrix<FieldData> >(cerr," \n")); cerr << endl;
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
	ublas::vector<FieldData,ublas::bounded_array<double,3> > normal(3);
	ierr = ShapeFaceNormalMBTRI(diffNTRI,coords_face,&*normal.data().begin()); CHKERRQ(ierr);
	double area = cblas_dnrm2(3,&*normal.data().begin(),1)*0.5;

	//higher order face shape
	vector< ublas::vector<FieldData> > Normals_at_Gauss_pts;
	ierr = GetHierarchicalGeometryApproximation_FaceNormal(siit->ent,Normals_at_Gauss_pts);  CHKERRQ(ierr);

	//calulate & assemble
	int g_dim = get_dim_gNTRI();
	for(int gg = 0;gg<g_dim;gg++) {
	  ublas::vector<FieldData,ublas::bounded_array<double,3> > traction_at_Gauss_pt = traction;
	  double w,area_at_Gauss_pt;
	  if(Normals_at_Gauss_pts.size()>0) {
	    if(Normals_at_Gauss_pts.size()!=(unsigned int)g_dim) {
	      SETERRQ(PETSC_COMM_SELF,1,"data inconsistency");
	    }
	    area_at_Gauss_pt = cblas_dnrm2(3,Normals_at_Gauss_pts[gg].data().begin(),1)*0.5;
	    w = area_at_Gauss_pt*G_W_TRI[gg];
	    traction_at_Gauss_pt += (pressure/(2*area_at_Gauss_pt))*Normals_at_Gauss_pts[gg];
	  } else {
	    w = area*G_W_TRI[gg];
	    traction_at_Gauss_pt += (pressure/area)*normal;
	  }

	  //nodes
	  ublas::vector<FieldData> f_ext_nodes(FaceNMatrix_nodes[0].size2());
	  f_ext_nodes = w*prod(trans(FaceNMatrix_nodes[gg]), traction_at_Gauss_pt);
	  ierr = VecSetValues(F_ext,RowGlob_nodes.size(),&(RowGlob_nodes)[0],&(f_ext_nodes.data())[0],ADD_VALUES); CHKERRQ(ierr);

	  //edges
	  siiit = side_table.get<1>().lower_bound(boost::make_tuple(MBEDGE,0));
	  for(;siiit!=hi_siiit;siiit++) {
	    vector<DofIdx>& RowGlob_edge = RowGlob[1+siiit->side_number];
	    if(RowGlob_edge.size()>0) {
	      vector<ublas::matrix<FieldData> >& FaceNMatrix_edge = FaceNMatrix_edges[siiit->side_number];
	      ublas::vector<FieldData> f_ext_edges = ublas::zero_vector<FieldData>(FaceNMatrix_edge[0].size2());
	      f_ext_edges = w*prod(trans(FaceNMatrix_edge[gg]), traction_at_Gauss_pt);
	      if(RowGlob_edge.size()!=f_ext_edges.size()) {
		SETERRQ(PETSC_COMM_SELF,1,"wrong size: RowGlob_edge.size()!=f_ext_edges.size()");
	      }
	      ierr = VecSetValues(F_ext,RowGlob_edge.size(),&(RowGlob_edge[0]),&(f_ext_edges.data())[0],ADD_VALUES); CHKERRQ(ierr);
	    }
	  }

	  //face
	  if(RowGlob_face.size()>0) {
	    ublas::vector<FieldData> f_ext_face(FaceNMatrix_face[0].size2());
	    f_ext_face = w*prod(trans(FaceNMatrix_face[gg]),traction_at_Gauss_pt);
	    if(RowGlob_face.size()!=f_ext_face.size()) SETERRQ(PETSC_COMM_SELF,1,"wrong size: RowGlob_face.size()!=f_ext_face.size()");
	    ierr = VecSetValues(F_ext,RowGlob_face.size(),&(RowGlob_face)[0],&(f_ext_face.data())[0],ADD_VALUES); CHKERRQ(ierr);
	  }
	}
      }	

      PetscFunctionReturn(0);

    }

    virtual PetscErrorCode NeumannBC(Vec F) {
      PetscFunctionBegin;
      
      for(_IT_CUBITMESHSETS_BY_BCDATA_TYPE_FOR_LOOP_(mField,NodeSet|ForceSet,it)) {
	
	ublas::vector<FieldData,ublas::bounded_array<double,3> > traction_glob(3,0);

	force_cubit_bc_data mydata;
	ierr = it->get_cubit_bc_data_structure(mydata); CHKERRQ(ierr);
	Range faces;
	ierr = it->get_Cubit_msId_entities_by_dimension(mField.get_moab(),2,faces,true); CHKERRQ(ierr);
	/*ostringstream ss;
	ss << *it << endl;
	ss << mydata;
	ss << "nb faces " << faces.size() << endl;
	PetscPrintf(PETSC_COMM_WORLD,ss.str().c_str());*/

	traction_glob[0] = mydata.data.value3;
	traction_glob[1] = mydata.data.value4;
	traction_glob[2] = mydata.data.value5;
	traction_glob *= mydata.data.value1;
  
	ierr = NeumannBC_Faces(F,0,traction_glob,faces); CHKERRQ(ierr);

      }

      ublas::vector<FieldData,ublas::bounded_array<double,3> > traction_glob = ublas::zero_vector<FieldData>(3);
      for(_IT_CUBITMESHSETS_BY_BCDATA_TYPE_FOR_LOOP_(mField,SideSet|PressureSet,it)) {

	pressure_cubit_bc_data mydata;
	ierr = it->get_cubit_bc_data_structure(mydata); CHKERRQ(ierr);
	Range faces;
	ierr = it->get_Cubit_msId_entities_by_dimension(mField.get_moab(),2,faces,true); CHKERRQ(ierr);

	ierr = NeumannBC_Faces(F,mydata.data.value1,traction_glob,faces); CHKERRQ(ierr);

      }

      PetscFunctionReturn(0);
    }

    virtual PetscErrorCode GetMatricesRows() {
      PetscFunctionBegin;
      //indicies ROWS
      row_mat = 0;
      ierr = GetRowGlobalIndices("DISPLACEMENT",RowGlob[row_mat]); CHKERRQ(ierr);
      ierr = GetRowLocalIndices("DISPLACEMENT",RowLocal[row_mat]); CHKERRQ(ierr);
      ierr = GetGaussRowNMatrix("DISPLACEMENT",rowNMatrices[row_mat]); CHKERRQ(ierr);
      ierr = GetGaussRowDiffNMatrix("DISPLACEMENT",rowDiffNMatrices[row_mat]); CHKERRQ(ierr);
      //HO gemometry
      ierr = GetHierarchicalGeometryApproximation_ApplyToDiffShapeFunction(3,invH,rowDiffNMatrices[row_mat]); CHKERRQ(ierr);
      //
      ierr = MakeBMatrix3D("DISPLACEMENT",rowDiffNMatrices[row_mat],rowBMatrices[row_mat]);  CHKERRQ(ierr);
      row_mat++;
      for(int ee = 0;ee<6;ee++) { //edges matrices
	ierr = GetRowGlobalIndices("DISPLACEMENT",MBEDGE,RowGlob[row_mat],ee); CHKERRQ(ierr);
	ierr = GetRowLocalIndices("DISPLACEMENT",MBEDGE,RowLocal[row_mat],ee); CHKERRQ(ierr);
	if(RowGlob[row_mat].size()!=0) {
	  ierr = GetGaussRowNMatrix("DISPLACEMENT",MBEDGE,rowNMatrices[row_mat],ee); CHKERRQ(ierr);
	  ierr = GetGaussRowDiffNMatrix("DISPLACEMENT",MBEDGE,rowDiffNMatrices[row_mat],ee); CHKERRQ(ierr);
	  //HO gemometry
	  ierr = GetHierarchicalGeometryApproximation_ApplyToDiffShapeFunction(3,invH,rowDiffNMatrices[row_mat]); CHKERRQ(ierr);
	  //
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
	  //HO gemometry
	  ierr = GetHierarchicalGeometryApproximation_ApplyToDiffShapeFunction(3,invH,rowDiffNMatrices[row_mat]); CHKERRQ(ierr);
	  //
	  ierr = MakeBMatrix3D("DISPLACEMENT",rowDiffNMatrices[row_mat],rowBMatrices[row_mat]);  CHKERRQ(ierr);
	  row_mat++;
	}
      }
      ierr = GetRowGlobalIndices("DISPLACEMENT",MBTET,RowGlob[row_mat]); CHKERRQ(ierr);
      ierr = GetRowLocalIndices("DISPLACEMENT",MBTET,RowLocal[row_mat]); CHKERRQ(ierr);
      if(RowGlob[row_mat].size() != 0) { //volume matrices
	ierr = GetGaussRowNMatrix("DISPLACEMENT",MBTET,rowNMatrices[row_mat]); CHKERRQ(ierr);
	ierr = GetGaussRowDiffNMatrix("DISPLACEMENT",MBTET,rowDiffNMatrices[row_mat]); CHKERRQ(ierr);
	//HO gemometry
	ierr = GetHierarchicalGeometryApproximation_ApplyToDiffShapeFunction(3,invH,rowDiffNMatrices[row_mat]); CHKERRQ(ierr);
	//
	ierr = MakeBMatrix3D("DISPLACEMENT",rowDiffNMatrices[row_mat],rowBMatrices[row_mat]);  CHKERRQ(ierr);
	row_mat++;
      }

      PetscFunctionReturn(0);
    }

    virtual PetscErrorCode GetMatricesCols() {
      PetscFunctionBegin;
      //Higher order approximation of geometry
      //ierr = GetHierarchicalGeometryApproximation(invH,detH); CHKERRQ(ierr);
      //indicies COLS
      col_mat = 0;
      ierr = GetColGlobalIndices("DISPLACEMENT",ColGlob[col_mat]); CHKERRQ(ierr);
      ierr = GetColLocalIndices("DISPLACEMENT",ColLocal[col_mat]); CHKERRQ(ierr);
      ierr = GetGaussColNMatrix("DISPLACEMENT",colNMatrices[col_mat]); CHKERRQ(ierr);
      ierr = GetGaussColDiffNMatrix("DISPLACEMENT",colDiffNMatrices[col_mat]); CHKERRQ(ierr);
      //HO gemometry
      ierr = GetHierarchicalGeometryApproximation_ApplyToDiffShapeFunction(3,invH,colDiffNMatrices[col_mat]); CHKERRQ(ierr);
      //
      ierr = MakeBMatrix3D("DISPLACEMENT",colDiffNMatrices[col_mat],colBMatrices[col_mat]);  CHKERRQ(ierr);
      col_mat++;
      for(int ee = 0;ee<6;ee++) { //edges matrices
	ierr = GetColGlobalIndices("DISPLACEMENT",MBEDGE,ColGlob[col_mat],ee); CHKERRQ(ierr);
	ierr = GetColLocalIndices("DISPLACEMENT",MBEDGE,ColLocal[col_mat],ee); CHKERRQ(ierr);
	if(ColGlob[col_mat].size()!=0) {
	  ierr = GetGaussColNMatrix("DISPLACEMENT",MBEDGE,colNMatrices[col_mat],ee); CHKERRQ(ierr);
	  ierr = GetGaussColDiffNMatrix("DISPLACEMENT",MBEDGE,colDiffNMatrices[col_mat],ee); CHKERRQ(ierr);
	  //HO gemometry
	  ierr = GetHierarchicalGeometryApproximation_ApplyToDiffShapeFunction(3,invH,colDiffNMatrices[col_mat]); CHKERRQ(ierr);
	  //
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
	  //HO gemometry
	  ierr = GetHierarchicalGeometryApproximation_ApplyToDiffShapeFunction(3,invH,colDiffNMatrices[col_mat]); CHKERRQ(ierr);
	  //
	  ierr = MakeBMatrix3D("DISPLACEMENT",colDiffNMatrices[col_mat],colBMatrices[col_mat]);  CHKERRQ(ierr);
	  col_mat++;
	}
      }
      ierr = GetColGlobalIndices("DISPLACEMENT",MBTET,ColGlob[col_mat]); CHKERRQ(ierr);
      ierr = GetColLocalIndices("DISPLACEMENT",MBTET,ColLocal[col_mat]); CHKERRQ(ierr);
      if(ColGlob[col_mat].size() != 0) { //volume matrices
	ierr = GetGaussColNMatrix("DISPLACEMENT",MBTET,colNMatrices[col_mat]); CHKERRQ(ierr);
	ierr = GetGaussColDiffNMatrix("DISPLACEMENT",MBTET,colDiffNMatrices[col_mat]); CHKERRQ(ierr);
	//HO gemometry
	ierr = GetHierarchicalGeometryApproximation_ApplyToDiffShapeFunction(3,invH,colDiffNMatrices[col_mat]); CHKERRQ(ierr);
	//
	ierr = MakeBMatrix3D("DISPLACEMENT",colDiffNMatrices[col_mat],colBMatrices[col_mat]);  CHKERRQ(ierr);
	col_mat++;
      }

      PetscFunctionReturn(0);
    }
   

    virtual PetscErrorCode GetMatrices() {
      PetscFunctionBegin;
      //Higher order approximation of geometry
      ierr = GetHierarchicalGeometryApproximation(invH,detH); CHKERRQ(ierr);
      //
      ierr = GetMatricesRows(); CHKERRQ(ierr);
      ierr = GetMatricesCols(); CHKERRQ(ierr);
      PetscFunctionReturn(0);
    }

    ublas::matrix<ublas::matrix<FieldData> > K;
    ublas::matrix<FieldData> BD;
    virtual PetscErrorCode Stiffness() {
      PetscFunctionBegin;

      double _lambda,_mu;
      ierr = GetMatParameters(&_lambda,&_mu); CHKERRQ(ierr);
      ierr = calulateD(_lambda,_mu); CHKERRQ(ierr);

      K.resize(row_mat,col_mat);
      int g_dim = g_NTET.size()/4;
      for(int rr = 0;rr<row_mat;rr++) {
	for(int gg = 0;gg<g_dim;gg++) {
	  ublas::matrix<FieldData> &row_Mat = (rowBMatrices[rr])[gg];
	  double w = V*G_W_TET[gg];
	  if(detH.size()>0) {
	    w *= detH[gg];
	  }
	  BD.resize(6,row_Mat.size2());
	  //ublas::noalias(BD) = prod( w*D,row_Mat );
	  cblas_dsymm(CblasRowMajor,CblasLeft,CblasUpper,
	    BD.size1(),BD.size2(),
	    w,&*D.data().begin(),D.size2(),
	    &*row_Mat.data().begin(),row_Mat.size2(),
	    0.,&*BD.data().begin(),BD.size2());
	  for(int cc = rr;cc<col_mat;cc++) {
	    ublas::matrix<FieldData> &col_Mat = (colBMatrices[cc])[gg];
	    if(gg == 0) {
	      K(rr,cc).resize(BD.size2(),col_Mat.size2());
	      //ublas::noalias(K(rr,cc)) = prod(trans(BD) , col_Mat ); // int BT*D*B
	      cblas_dgemm(CblasRowMajor,CblasTrans,CblasNoTrans,
		BD.size2(),col_Mat.size2(),BD.size1(),
		1.,&*BD.data().begin(),BD.size2(),
		&*col_Mat.data().begin(),col_Mat.size2(),
		0.,&*K(rr,cc).data().begin(),K(rr,cc).size2());
	    } else {
	      //ublas::noalias(K(rr,cc)) += prod(trans(BTD) , col_Mat ); // int BT*D*B
	      cblas_dgemm(CblasRowMajor,CblasTrans,CblasNoTrans,
		BD.size2(),col_Mat.size2(),BD.size1(),
		1.,&*BD.data().begin(),BD.size2(),
		&*col_Mat.data().begin(),col_Mat.size2(),
		1.,&*K(rr,cc).data().begin(),K(rr,cc).size2());
	    }
	  }
	}
      }
      PetscFunctionReturn(0);
    }

    virtual PetscErrorCode Lhs() {
      PetscFunctionBegin;
      ierr = Stiffness(); CHKERRQ(ierr);
      for(int rr = 0;rr<row_mat;rr++) {
	for(int cc = rr;cc<col_mat;cc++) {
	  if(ColGlob[cc].size()==0) continue;
	  if(RowGlob[rr].size()!=K(rr,cc).size1()) SETERRQ(PETSC_COMM_SELF,1,"data inconsistency");
	  if(ColGlob[cc].size()!=K(rr,cc).size2()) SETERRQ(PETSC_COMM_SELF,1,"data inconsistency");
	  ierr = MatSetValues(Aij,RowGlob[rr].size(),&(RowGlob[rr])[0],ColGlob[cc].size(),&(ColGlob[cc])[0],&(K(rr,cc).data())[0],ADD_VALUES); CHKERRQ(ierr);
	  if(rr!=cc) {
	    K(cc,rr) = trans(K(rr,cc));
	    ierr = MatSetValues(Aij,ColGlob[cc].size(),&(ColGlob[cc])[0],RowGlob[rr].size(),&(RowGlob[rr])[0],&(K(cc,rr).data())[0],ADD_VALUES); CHKERRQ(ierr);
	  }
	}
      }
      PetscFunctionReturn(0);
    }

    virtual PetscErrorCode Rhs() {
      PetscFunctionBegin;
      ierr = Fint(F); CHKERRQ(ierr);
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

    vector<ublas::vector<FieldData> > f_int;
    virtual PetscErrorCode Fint() {
      PetscFunctionBegin;

      try {

      //Higher order approximation of geometry
      ierr = GetHierarchicalGeometryApproximation(invH,detH); CHKERRQ(ierr);

      double _lambda,_mu;
      ierr = GetMatParameters(&_lambda,&_mu); CHKERRQ(ierr);
      ierr = calulateD(_lambda,_mu); CHKERRQ(ierr);

      //Gradient at Gauss points; 
      vector< ublas::matrix< FieldData > > GradU_at_GaussPt;
      ierr = GetGaussDiffDataVector("DISPLACEMENT",GradU_at_GaussPt); CHKERRQ(ierr);
      unsigned int g_dim = g_NTET.size()/4;
      assert(GradU_at_GaussPt.size() == g_dim);
      NOT_USED(g_dim);
      vector< ublas::matrix< FieldData > >::iterator viit = GradU_at_GaussPt.begin();
      int gg = 0;
      for(;viit!=GradU_at_GaussPt.end();viit++,gg++) {
	try {
	  ublas::matrix< FieldData > GradU = *viit;
	  if(!invH.empty()) {
	    //GradU = 
	      //[ dU/dChi1 dU/dChi2 dU/dChi3 ]
	      //[ dV/dChi1 dV/dChi2 dU/dChi3 ]
	      //[ dW/dChi1 dW/dChi2 dW/dChi3 ]
	    //H = 
	      //[ dX1/dChi1 dX1/dChi2 dX1/dChi3 ]
	      //[ dX2/dChi1 dX2/dChi2 dX2/dChi3 ]
	      //[ dX3/dChi1 dX3/dChi2 dX3/dChi3 ]    
	    //invH = 
	      //[ dChi1/dX1 dChi1/dX2 dChi1/dX3 ]
	      //[ dChi2/dX1 dChi2/dX2 dChi2/dX3 ]
	      //[ dChi3/dX1 dChi3/dX2 dChi3/dX3 ]
	    //GradU = 
	      //[ dU/dX1 dU/dX2 dU/dX3 ]
	      //[ dV/dX1 dV/dX2 dV/dX3 ] = GradU * invH
	      //[ dW/dX1 dW/dX2 dW/dX3 ] 
	    GradU = prod( GradU, invH[gg] ); 
	  }
	  ublas::matrix< FieldData > Strain = 0.5*( GradU + trans(GradU) );
	  ublas::vector< FieldData > VoightStrain(6);
	  VoightStrain[0] = Strain(0,0);
	  VoightStrain[1] = Strain(1,1);
	  VoightStrain[2] = Strain(2,2);
	  VoightStrain[3] = 2*Strain(0,1);
	  VoightStrain[4] = 2*Strain(1,2);
	  VoightStrain[5] = 2*Strain(2,0);
	  double w = V*G_W_TET[gg];
	  ublas::vector<FieldData> VoightStress = prod(w*D,VoightStrain);
	  //BT * VoigtStress
	  for(int rr = 0;rr<row_mat;rr++) {
	    f_int.resize(row_mat);
	    ublas::matrix<FieldData> &B = (rowBMatrices[rr])[gg];
	    if(gg == 0) {
	      f_int[rr] = prod( trans(B), VoightStress );
	    } else {
	      f_int[rr] += prod( trans(B), VoightStress );
	    }
	  }
	} catch (const std::exception& ex) {
	  ostringstream ss;
	  ss << "thorw in method: " << ex.what() << " at line " << __LINE__ << " in file " << __FILE__ << endl;
	  SETERRQ(PETSC_COMM_SELF,1,ss.str().c_str());
	} 
      }

      } catch (const std::exception& ex) {
	ostringstream ss;
	ss << "thorw in method: " << ex.what() << " at line " << __LINE__ << " in file " << __FILE__ << endl;
	SETERRQ(PETSC_COMM_SELF,1,ss.str().c_str());
      }

      PetscFunctionReturn(0);
    }

    virtual PetscErrorCode Fint(Vec F_int) {
      PetscFunctionBegin;
      try {
	ierr = Fint(); CHKERRQ(ierr);
	for(int rr = 0;rr<row_mat;rr++) {
	  if(RowGlob[rr].size()!=f_int[rr].size()) SETERRQ(PETSC_COMM_SELF,1,"data inconsistency");
	  if(RowGlob[rr].size()==0) continue;
	  ierr = VecSetValues(F_int,RowGlob[rr].size(),&(RowGlob[rr])[0],&(f_int[rr].data()[0]),ADD_VALUES); CHKERRQ(ierr);
	}
      } catch (const std::exception& ex) {
	ostringstream ss;
	ss << "thorw in method: " << ex.what() << " at line " << __LINE__ << " in file " << __FILE__ << endl;
	SETERRQ(PETSC_COMM_SELF,1,ss.str().c_str());
      }
      PetscFunctionReturn(0);
    }


    virtual PetscErrorCode RhsAndLhs() {
      PetscFunctionBegin;

      ierr = Rhs(); CHKERRQ(ierr);
      ierr = Lhs(); CHKERRQ(ierr);

      PetscFunctionReturn(0);
    }

    PetscErrorCode operator()() {
      PetscFunctionBegin;
      ierr = OpStudentStart_TET(g_NTET); CHKERRQ(ierr);
      ierr = GetMatrices(); CHKERRQ(ierr);
      //Dirihlet Boundary Condition
      ierr = dirihlet_bc_method_ptr->SetDirihletBC_to_ElementIndicies(this,RowGlob,ColGlob,DirihletBC); CHKERRQ(ierr);

      //Assembly Aij and F
      ierr = RhsAndLhs(); CHKERRQ(ierr);

      //Neumann Boundary Conditions
      ierr = NeumannBC(F); CHKERRQ(ierr);

      ierr = OpStudentEnd(); CHKERRQ(ierr);
      PetscFunctionReturn(0);
    }

};

    
}

#endif //__ELASTICFEMETHOD_HPP__
