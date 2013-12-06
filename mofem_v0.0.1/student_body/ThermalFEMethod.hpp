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

#ifndef __THERMALFEMETHOD_HPP__
#define __THERMALFEMETHOD_HPP__

#include <boost/numeric/ublas/symmetric.hpp>

namespace MoFEM {

struct ThermalFEMethod: public FEMethod_UpLevelStudent {

    FieldInterface& mField;
    Mat Aij;
    Vec Data,F;

    ThermalFEMethod(
      FieldInterface& _mField): FEMethod_UpLevelStudent(_mField.get_moab(),1), mField(_mField),
      Aij(PETSC_NULL),Data(PETSC_NULL),F(PETSC_NULL) {};

    ThermalFEMethod(
      FieldInterface& _mField,BaseDirihletBC *_dirihlet_ptr,Mat &_Aij,Vec &_D,Vec& _F,
      double _Ther_Cond):
      FEMethod_UpLevelStudent(_mField.get_moab(),_dirihlet_ptr,1), mField(_mField),
      Aij(_Aij),Data(_D),F(_F),
                    Ther_Cond(_Ther_Cond){

      pcomm = ParallelComm::get_pcomm(&moab,MYPCOMM_INDEX);
        
      RowGlob.resize(1+6+4+1);
      RowLocal.resize(1+6+4+1);
      rowNMatrices.resize(1+6+4+1);
      rowDiffNMatrices.resize(1+6+4+1);
      //rowBMatrices.resize(1+6+4+1);
      ColGlob.resize(1+6+4+1);
      ColLocal.resize(1+6+4+1);
      colNMatrices.resize(1+6+4+1);
      colDiffNMatrices.resize(1+6+4+1);
      //colBMatrices.resize(1+6+4+1);

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
      
      
    }; 

    ErrorCode rval;
    
    ParallelComm* pcomm;
    PetscLogDouble t1,t2;
    PetscLogDouble v1,v2;

    double Ther_Cond;
    ublas::matrix<FieldData> D_lambda,D_mu,D;

    int row_mat,col_mat;
    vector<vector<DofIdx> > RowGlob;
    vector<vector<DofIdx> > RowLocal;
    vector<vector<ublas::matrix<FieldData> > > rowNMatrices;
    vector<vector<ublas::matrix<FieldData> > > rowDiffNMatrices;
    //vector<vector<ublas::matrix<FieldData> > > rowBMatrices;
    vector<vector<DofIdx> > ColGlob;
    vector<vector<DofIdx> > ColLocal;
    vector<vector<ublas::matrix<FieldData> > > colNMatrices;
    vector<vector<ublas::matrix<FieldData> > > colDiffNMatrices;
    //vector<vector<ublas::matrix<FieldData> > > colBMatrices;

    vector<DofIdx> DirihletBC;

    vector<double> g_NTET,g_NTRI;
    const double* G_W_TET;
    const double* G_W_TRI;

    virtual PetscErrorCode calulateD(double _Ther_Cond) {
      PetscFunctionBegin;

        //cout<<"D_lambda "<<D_lambda<<endl;
      D = _Ther_Cond*D_lambda;
        //cout<<"D "<<D<<endl;
        
      PetscFunctionReturn(0);
    }
    
    PetscErrorCode GetMatParameters(double *_Ther_Cond) {
        PetscFunctionBegin;
        
        *_Ther_Cond = Ther_Cond;
        
        EntityHandle ent = fe_ptr->get_ent();
        for(_IT_CUBITMESHSETS_BY_BCDATA_TYPE_FOR_LOOP_(mField,BlockSet|Mat_ThermalSet,it)) {
            
            Mat_Thermal mydata;
            ierr = it->get_attribute_data_structure(mydata); CHKERRQ(ierr);
            
            Range meshsets;
            rval = moab.get_entities_by_type(it->meshset,MBENTITYSET,meshsets,true); CHKERR_PETSC(rval);
            for(Range::iterator mit = meshsets.begin();mit != meshsets.end(); mit++) {
                if( moab.contains_entities(*mit,&ent,1) ) {
                    *_Ther_Cond = mydata.data.Conductivity;
                    //cout<< " mydata.data.Conductivity "<<mydata.data.Conductivity<<endl;
                    break;
                }
            }
        }
        PetscFunctionReturn(0);
    }

    PetscErrorCode preProcess() {
      PetscFunctionBegin;
      PetscSynchronizedPrintf(PETSC_COMM_WORLD,"Start Assembly\n");
      ierr = PetscTime(&v1); CHKERRQ(ierr);
      ierr = PetscGetCPUTime(&t1); CHKERRQ(ierr);
      g_NTET.resize(4*45);
      ShapeMBTET(&g_NTET[0],G_TET_X45,G_TET_Y45,G_TET_Z45,45);
      G_W_TET = G_TET_W45;
      g_NTRI.resize(3*13);
      ShapeMBTRI(&g_NTRI[0],G_TRI_X13,G_TRI_Y13,13); 
      G_W_TRI = G_TRI_W13;

      // See FEAP - - A Finite Element Analysis Program
      D_lambda.resize(3,3);
      D_lambda.clear();
      for(int rr = 0;rr<3;rr++) {
	for(int cc = 0;cc<3;cc++) {
        if(rr==cc) D_lambda(rr,cc) = 1;
	}
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
      PetscFunctionReturn(0);
    }

    virtual PetscErrorCode NeumannBC_Faces(Vec F_ext,ublas::vector<FieldData,ublas::bounded_array<double,1> > &flux,Range& faces) {
      PetscFunctionBegin;

      SideNumber_multiIndex& side_table = const_cast<SideNumber_multiIndex&>(fe_ent_ptr->get_side_number_table());
      SideNumber_multiIndex::nth_index<1>::type::iterator siit = side_table.get<1>().lower_bound(boost::make_tuple(MBTRI,0));
      SideNumber_multiIndex::nth_index<1>::type::iterator hi_siit = side_table.get<1>().upper_bound(boost::make_tuple(MBTRI,4));
      
        //cout<<"NeumannBC "<<endl;

      for(;siit!=hi_siit;siit++) {
	Range::iterator fit = find(faces.begin(),faces.end(),siit->ent);
	if(fit==faces.end()) continue;

	ierr = ShapeFunctions_TRI(siit->ent,g_NTRI);  CHKERRQ(ierr);

	//nodes
	vector<DofIdx>& RowGlob_nodes = RowGlob[0];
	vector< ublas::matrix<FieldData> > FaceNMatrix_nodes;
	ierr = GetGaussRowFaceNMatrix(siit->ent,"TEMPERATURE",FaceNMatrix_nodes,MBVERTEX); CHKERRQ(ierr);
    //cout<<"FaceNMatrix_nodes  "<<FaceNMatrix_nodes[0]<<endl;
 	//copy(FaceNMatrix_nodes.begin(),FaceNMatrix_nodes.end(),ostream_iterator<ublas::matrix<FieldData> >(cerr," \n")); cerr << endl;

	//edges
	vector<vector<ublas::matrix<FieldData> > > FaceNMatrix_edges(6);
	SideNumber_multiIndex::nth_index<1>::type::iterator siiit = side_table.get<1>().lower_bound(boost::make_tuple(MBEDGE,0));
	SideNumber_multiIndex::nth_index<1>::type::iterator hi_siiit = side_table.get<1>().upper_bound(boost::make_tuple(MBEDGE,6));
	for(;siiit!=hi_siiit;siiit++) {
	  if(RowGlob[1+siiit->side_number].size()>0) {
	    ierr = GetGaussRowFaceNMatrix(siit->ent,"TEMPERATURE",FaceNMatrix_edges[siiit->side_number],MBEDGE,siiit->ent); CHKERRQ(ierr);
	    //cerr << "ee ";
	    //copy(FaceNMatrix_face.begin(),FaceNMatrix_face.end(),ostream_iterator<ublas::matrix<FieldData> >(cerr," \n")); cerr << endl;
	  }
	}

	//faces
	vector<DofIdx>& RowGlob_face = RowGlob[1+6+siit->side_number];
	vector< ublas::matrix<FieldData> > FaceNMatrix_face;
	if(RowGlob_face.size()>0) {
	  ierr = GetGaussRowFaceNMatrix(siit->ent,"TEMPERATURE",FaceNMatrix_face,MBTRI); CHKERRQ(ierr);
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
    //cout<<"area "<<area<<endl;

	//calulate & assemble

	//nodes
	ublas::vector<FieldData> f_ext_nodes = ublas::zero_vector<FieldData>(FaceNMatrix_nodes[0].size2());
	int g_dim = get_dim_gNTRI();
	for(int gg = 0;gg<g_dim;gg++) {
	  double w = area*G_W_TRI[gg];
	  f_ext_nodes += w*prod(trans(FaceNMatrix_nodes[gg]), flux);
	}
    //cout<<"f_ext_nodes.size()  "<<f_ext_nodes<<endl;
    //cout<<"RowGlob_nodes.size()  "<<RowGlob_nodes.size()<<endl;
    //cout<<"RowGlob_nodes "<< RowGlob_nodes<<endl;
    //cout<<"f_ext_nodes.data() "<< f_ext_nodes.data()<<endl;
      
          
	ierr = VecSetValues(F_ext,RowGlob_nodes.size(),&(RowGlob_nodes)[0],&(f_ext_nodes.data())[0],ADD_VALUES); CHKERRQ(ierr);
    //cout<<"F_ext  "<<F_ext.size()<<endl;
    //cout<<"F_ext  "<<F_ext <<endl;
    
	//edges
	siiit = side_table.get<1>().lower_bound(boost::make_tuple(MBEDGE,0));
	for(;siiit!=hi_siiit;siiit++) {
	  vector<DofIdx>& RowGlob_edge = RowGlob[1+siiit->side_number];
	  if(RowGlob_edge.size()>0) {
	    vector<ublas::matrix<FieldData> >& FaceNMatrix_edge = FaceNMatrix_edges[siiit->side_number];
	    ublas::vector<FieldData> f_ext_edges = ublas::zero_vector<FieldData>(FaceNMatrix_edge[0].size2());
	    for(int gg = 0;gg<g_dim;gg++) {
	      double w = area*G_W_TRI[gg];
	      f_ext_edges += w*prod(trans(FaceNMatrix_edge[gg]), flux);
	    }
	    if(RowGlob_edge.size()!=f_ext_edges.size()) SETERRQ(PETSC_COMM_SELF,1,"wrong size: RowGlob_edge.size()!=f_ext_edges.size()");
	    ierr = VecSetValues(F_ext,RowGlob_edge.size(),&(RowGlob_edge[0]),&(f_ext_edges.data())[0],ADD_VALUES); CHKERRQ(ierr);
	  }
	}
  
 	//face     
	if(RowGlob_face.size()>0) {
	  ublas::vector<FieldData> f_ext_faces = ublas::zero_vector<FieldData>(FaceNMatrix_face[0].size2());
	  for(int gg = 0;gg<g_dim;gg++) {
	    double w = area*G_W_TRI[gg];
	    f_ext_faces += w*prod(trans(FaceNMatrix_face[gg]), flux);
	  }
	  if(RowGlob_face.size()!=f_ext_faces.size()) SETERRQ(PETSC_COMM_SELF,1,"wrong size: RowGlob_face.size()!=f_ext_faces.size()");
	  ierr = VecSetValues(F_ext,RowGlob_face.size(),&(RowGlob_face)[0],&(f_ext_faces.data())[0],ADD_VALUES); CHKERRQ(ierr);
	}

      }

      PetscFunctionReturn(0);

    }

    virtual PetscErrorCode NeumannBC(Vec F) {
      PetscFunctionBegin;
      
      for(_IT_CUBITMESHSETS_BY_BCDATA_TYPE_FOR_LOOP_(mField,SideSet|HeatfluxSet,it)) {
	
	ublas::vector<FieldData,ublas::bounded_array<double,1> > flux(1);

	heatflux_cubit_bc_data mydata;
	ierr = it->get_cubit_bc_data_structure(mydata); CHKERRQ(ierr);
	Range faces;
	ierr = it->get_Cubit_msId_entities_by_dimension(mField.get_moab(),2,faces,true); CHKERRQ(ierr);
	/*ostringstream ss;
	ss << *it << endl;
	ss << mydata;
	ss << "nb faces " << faces.size() << endl;
	PetscPrintf(PETSC_COMM_WORLD,ss.str().c_str());*/

	flux[0] = mydata.data.value1;
	//cout<<"flux "<<flux[0]<<endl;

	ierr = NeumannBC_Faces(F,flux,faces); CHKERRQ(ierr);

      }

      PetscFunctionReturn(0);
    }

    virtual PetscErrorCode GetMatricesRows() {
      PetscFunctionBegin;
      //indicies ROWS
      //cout<<"GetMatricesRows " <<endl;
      row_mat = 0;
      ierr = GetRowGlobalIndices("TEMPERATURE",RowGlob[row_mat]); CHKERRQ(ierr);    //local belong to processors
      ierr = GetRowLocalIndices("TEMPERATURE",RowLocal[row_mat]); CHKERRQ(ierr);
      ierr = GetGaussRowNMatrix("TEMPERATURE",rowNMatrices[row_mat]); CHKERRQ(ierr);
      ierr = GetGaussRowDiffNMatrix("TEMPERATURE",rowDiffNMatrices[row_mat]); CHKERRQ(ierr);
        
        //cout<<"RowLocal "<< (RowLocal[0])[0]<<"\t"<<(RowLocal[0])[1]<<"\t"<<(RowLocal[0])[2]<<"\t"<<(RowLocal[0])[3]<<"\t"<<endl;
        //cout<<"RowGlob  "<< (RowGlob[0])[0]<<"\t"<<(RowGlob[0])[1]<<"\t"<<(RowGlob[0])[2]<<"\t"<<(RowGlob[0])[3]<<"\t"<<(RowGlob[0])[4]<<"\t"<<(RowGlob[0])[5]<<"\t"<<endl;

     //ierr = MakeBMatrix3D("TEMPERATURE",rowDiffNMatrices[row_mat],rowBMatrices[row_mat]);  CHKERRQ(ierr);
        
        
      row_mat++;
      for(int ee = 0;ee<6;ee++) { //edges matrices
	ierr = GetRowGlobalIndices("TEMPERATURE",MBEDGE,RowGlob[row_mat],ee); CHKERRQ(ierr);
	ierr = GetRowLocalIndices("TEMPERATURE",MBEDGE,RowLocal[row_mat],ee); CHKERRQ(ierr);
	if(RowGlob[row_mat].size()!=0) {
	  ierr = GetGaussRowNMatrix("TEMPERATURE",MBEDGE,rowNMatrices[row_mat],ee); CHKERRQ(ierr);
	  ierr = GetGaussRowDiffNMatrix("TEMPERATURE",MBEDGE,rowDiffNMatrices[row_mat],ee); CHKERRQ(ierr);
	  //ierr = MakeBMatrix3D("TEMPERATURE",rowDiffNMatrices[row_mat],rowBMatrices[row_mat]);  CHKERRQ(ierr);
	  row_mat++;
	}
      }
      for(int ff = 0;ff<4;ff++) { //faces matrices
	ierr = GetRowGlobalIndices("TEMPERATURE",MBTRI,RowGlob[row_mat],ff); CHKERRQ(ierr);
	ierr = GetRowLocalIndices("TEMPERATURE",MBTRI,RowLocal[row_mat],ff); CHKERRQ(ierr);
	if(RowGlob[row_mat].size()!=0) {
	  ierr = GetGaussRowNMatrix("TEMPERATURE",MBTRI,rowNMatrices[row_mat],ff); CHKERRQ(ierr);
	  ierr = GetGaussRowDiffNMatrix("TEMPERATURE",MBTRI,rowDiffNMatrices[row_mat],ff); CHKERRQ(ierr);
	  //ierr = MakeBMatrix3D("TEMPERATURE",rowDiffNMatrices[row_mat],rowBMatrices[row_mat]);  CHKERRQ(ierr);
	  row_mat++;
	}
      }
      ierr = GetRowGlobalIndices("TEMPERATURE",MBTET,RowGlob[row_mat]); CHKERRQ(ierr);
      ierr = GetRowLocalIndices("TEMPERATURE",MBTET,RowLocal[row_mat]); CHKERRQ(ierr);
      if(RowGlob[row_mat].size() != 0) { //volume matrices
	ierr = GetGaussRowNMatrix("TEMPERATURE",MBTET,rowNMatrices[row_mat]); CHKERRQ(ierr);
	ierr = GetGaussRowDiffNMatrix("TEMPERATURE",MBTET,rowDiffNMatrices[row_mat]); CHKERRQ(ierr);
	//ierr = MakeBMatrix3D("TEMPERATURE",rowDiffNMatrices[row_mat],rowBMatrices[row_mat]);  CHKERRQ(ierr);
	row_mat++;
      }
      PetscFunctionReturn(0);
    }

    virtual PetscErrorCode GetMatricesCols() {
      PetscFunctionBegin;
      //indicies COLS
      //cout<<"GetMatricesCols " <<endl;
      col_mat = 0;
      ierr = GetColGlobalIndices("TEMPERATURE",ColGlob[col_mat]); CHKERRQ(ierr);
      ierr = GetColLocalIndices("TEMPERATURE",ColLocal[col_mat]); CHKERRQ(ierr);
      ierr = GetGaussColNMatrix("TEMPERATURE",colNMatrices[col_mat]); CHKERRQ(ierr);
      ierr = GetGaussColDiffNMatrix("TEMPERATURE",colDiffNMatrices[col_mat]); CHKERRQ(ierr);
      //cout<<"colDiffNMatrices "<< (colDiffNMatrices[0])[0] <<endl;
      //ierr = MakeBMatrix3D("TEMPERATURE",colDiffNMatrices[col_mat],colBMatrices[col_mat]);  CHKERRQ(ierr);
      col_mat++;
      for(int ee = 0;ee<6;ee++) { //edges matrices
	ierr = GetColGlobalIndices("TEMPERATURE",MBEDGE,ColGlob[col_mat],ee); CHKERRQ(ierr);
	ierr = GetColLocalIndices("TEMPERATURE",MBEDGE,ColLocal[col_mat],ee); CHKERRQ(ierr);
	if(ColGlob[col_mat].size()!=0) {
	  ierr = GetGaussColNMatrix("TEMPERATURE",MBEDGE,colNMatrices[col_mat],ee); CHKERRQ(ierr);
	  ierr = GetGaussColDiffNMatrix("TEMPERATURE",MBEDGE,colDiffNMatrices[col_mat],ee); CHKERRQ(ierr);
      //cout<<"colDiffNMatrices Edges "<< (colDiffNMatrices[0])[0] <<endl;
	  //ierr = MakeBMatrix3D("TEMPERATURE",colDiffNMatrices[col_mat],colBMatrices[col_mat]);  CHKERRQ(ierr);
	  col_mat++;
	}
      }
      for(int ff = 0;ff<4;ff++) { //faces matrices
	ierr = GetColGlobalIndices("TEMPERATURE",MBTRI,ColGlob[col_mat],ff); CHKERRQ(ierr);
	ierr = GetColLocalIndices("TEMPERATURE",MBTRI,ColLocal[col_mat],ff); CHKERRQ(ierr);
	if(ColGlob[col_mat].size()!=0) {
	  ierr = GetGaussColNMatrix("TEMPERATURE",MBTRI,colNMatrices[col_mat],ff); CHKERRQ(ierr);
	  ierr = GetGaussColDiffNMatrix("TEMPERATURE",MBTRI,colDiffNMatrices[col_mat],ff); CHKERRQ(ierr);
	  //ierr = MakeBMatrix3D("TEMPERATURE",colDiffNMatrices[col_mat],colBMatrices[col_mat]);  CHKERRQ(ierr);
	  col_mat++;
	}
      }
      ierr = GetColGlobalIndices("TEMPERATURE",MBTET,ColGlob[col_mat]); CHKERRQ(ierr);
      ierr = GetColLocalIndices("TEMPERATURE",MBTET,ColLocal[col_mat]); CHKERRQ(ierr);
      if(ColGlob[col_mat].size() != 0) { //volume matrices
	ierr = GetGaussColNMatrix("TEMPERATURE",MBTET,colNMatrices[col_mat]); CHKERRQ(ierr);
	ierr = GetGaussColDiffNMatrix("TEMPERATURE",MBTET,colDiffNMatrices[col_mat]); CHKERRQ(ierr);
	//ierr = MakeBMatrix3D("TEMPERATURE",colDiffNMatrices[col_mat],colBMatrices[col_mat]);  CHKERRQ(ierr);
	col_mat++;
      }
      PetscFunctionReturn(0);
    }
   
    virtual PetscErrorCode GetMatrices() {
      PetscFunctionBegin;
      //cout<< "GetMatrices " << endl;
      ierr = GetMatricesRows(); CHKERRQ(ierr);
      ierr = GetMatricesCols(); CHKERRQ(ierr);
      PetscFunctionReturn(0);
    }

    ublas::matrix<ublas::matrix<FieldData> > K;
    ublas::matrix<FieldData> BD;
    virtual PetscErrorCode Stiffness() {
      PetscFunctionBegin;
      //cout<<"Stiffness "<<endl;

      double Ther_Cond;
      ierr = GetMatParameters(&Ther_Cond); CHKERRQ(ierr);
      //cout<< "Ther_Cond "<<Ther_Cond<<endl;
      ierr = calulateD(Ther_Cond); CHKERRQ(ierr);

      //cout<< "D "<<D<<endl;
      K.resize(row_mat,col_mat);
      int g_dim = g_NTET.size()/4;
        for(int rr = 0;rr<row_mat;rr++) {
	for(int gg = 0;gg<g_dim;gg++) {
	  ublas::matrix<FieldData> &row_Mat = (rowDiffNMatrices[rr])[gg];
	  double w = V*G_W_TET[gg];
        
      //cout<<"row_Mat.size2() "<<(rowDiffNMatrices[rr])[gg]<<endl;
        //cout<<"D "<<D<<endl;
        //cout<<"w "<<w<<endl;
	  BD.resize(3,row_Mat.size2());
	  //ublas::noalias(BD) = prod( w*D,row_Mat );
	  cblas_dsymm(CblasRowMajor,CblasLeft,CblasUpper,
	    BD.size1(),BD.size2(),
	    w,&*D.data().begin(),D.size2(),
	    &*row_Mat.data().begin(),row_Mat.size2(),
	    0.,&*BD.data().begin(),BD.size2());
        
        //cout<<"BD "<<BD<<endl;
        
	  for(int cc = rr;cc<col_mat;cc++) {
	    ublas::matrix<FieldData> &col_Mat = (rowDiffNMatrices[cc])[gg];
	    if(gg == 0) {     //for first gauss point k= BT*D*B + 0*k
	      K(rr,cc).resize(BD.size2(),col_Mat.size2());
	      //ublas::noalias(K(rr,cc)) = prod(trans(BD) , col_Mat ); // int BT*D*B
	      cblas_dgemm(CblasRowMajor,CblasTrans,CblasNoTrans,
		BD.size2(),col_Mat.size2(),BD.size1(),
		1.,&*BD.data().begin(),BD.size2(),
		&*col_Mat.data().begin(),col_Mat.size2(),
		0.,&*K(rr,cc).data().begin(),K(rr,cc).size2());
	    } else {   // for gauss other than first  k= BT*D*B + 1*k
	      //ublas::noalias(K(rr,cc)) += prod(trans(BTD) , col_Mat ); // int BT*D*B
	      cblas_dgemm(CblasRowMajor,CblasTrans,CblasNoTrans,
		BD.size2(),col_Mat.size2(),BD.size1(),
		1.,&*BD.data().begin(),BD.size2(),
		&*col_Mat.data().begin(),col_Mat.size2(),
		1.,&*K(rr,cc).data().begin(),K(rr,cc).size2());
            //cout<< "K  "<<K<<endl;
	    }
	  }
	}
     }
        //cout<< "K1  "<<K<<endl;
        //cout<<endl<<endl;

        
       PetscFunctionReturn(0);
    }

    virtual PetscErrorCode Lhs() {
      PetscFunctionBegin;
      //cout<<"Lhs "<<endl;
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
        
    
//        cout<<"F1 "<<endl<<endl<<endl<<endl;
//        PetscViewer    viewer;
//        VecView(F,viewer);
//        cout<<endl<<endl<<endl<<endl;
        
        
       ierr = Fint(F); CHKERRQ(ierr);
        
//        cout<<"F2 "<<endl<<endl<<endl<<endl;
//        VecView(F,viewer);
//        cout<<endl<<endl<<endl<<endl;

        
        
        
      ublas::vector<FieldData> f[row_mat];
      int g_dim = g_NTET.size()/4;
        for(int rr = 0;rr<row_mat;rr++) {
	for(int gg = 0;gg<g_dim;gg++) {
	  ublas::matrix<FieldData> &row_Mat = (rowNMatrices[rr])[gg];
        if(gg == 0) f[rr] = ublas::zero_vector<FieldData>(row_Mat.size2());  //cout<<"row_Mat.size2()"<<row_Mat.size2()<<endl;
	  ///f matrices
	  // calulate body force (f.e. garvity force
	}
	if(RowGlob[rr].size()!=f[rr].size()) SETERRQ(PETSC_COMM_SELF,1,"data inconsistency");
	if(RowGlob[rr].size()==0) continue;
	ierr = VecSetValues(F,RowGlob[rr].size(),&(RowGlob[rr])[0],&(f[rr].data())[0],ADD_VALUES); CHKERRQ(ierr);
      }
        
        //cout<<"Rhs "<<endl<<endl<<endl<<endl;
        //PetscViewer    viewer;
        //VecView(F,viewer);
        //cout<<endl<<endl<<endl<<endl;
        

      PetscFunctionReturn(0);
        
    }

    vector<ublas::vector<FieldData> > f_int;
    virtual PetscErrorCode Fint() {
      PetscFunctionBegin;

      double Ther_Cond;
      ierr = GetMatParameters(&Ther_Cond); CHKERRQ(ierr);
      ierr = calulateD(Ther_Cond); CHKERRQ(ierr);
        //cout<< "Ther_Cond "<<Ther_Cond<<endl;
        //cout<< "D "<<D<<endl;
      //Gradient at Gauss points;
      vector< ublas::matrix< FieldData > > GradU_at_GaussPt;
      ierr = GetGaussDiffDataVector("TEMPERATURE",GradU_at_GaussPt); CHKERRQ(ierr);
       
        
      unsigned int g_dim = g_NTET.size()/4;
      assert(GradU_at_GaussPt.size() == g_dim);
      NOT_USED(g_dim);
      vector< ublas::matrix< FieldData > >::iterator viit = GradU_at_GaussPt.begin();
      int gg = 0;
      for(;viit!=GradU_at_GaussPt.end();viit++,gg++) {
	  ublas::matrix< FieldData > GradU = *viit;
      //cout<<"GradU "<< GradU << endl;
	  ublas::matrix< FieldData > Strain = 0.5*( GradU + GradU );
      //cout<<"Strain "<< Strain << endl;
	  ublas::vector< FieldData > VoightStrain(3);
	  VoightStrain[0] = Strain(0,0);
	  VoightStrain[1] = Strain(0,1);
	  VoightStrain[2] = Strain(0,2);
      //cout<<"VoightStrain "<< VoightStrain << endl;
	  double w = V*G_W_TET[gg];
	  ublas::vector<FieldData> VoightStress = prod(w*D,VoightStrain);
	  //BT * VoigtStress
	  for(int rr = 0;rr<row_mat;rr++) {
	    f_int.resize(row_mat);
	    ublas::matrix<FieldData> &B = (rowDiffNMatrices[rr])[gg];
	    if(gg == 0) {
	      f_int[rr] = prod( trans(B), VoightStress );
	    } else {
	      f_int[rr] += prod( trans(B), VoightStress );
	    }
	  }
      }

      PetscFunctionReturn(0);
    }

    virtual PetscErrorCode Fint(Vec F_int) {
      PetscFunctionBegin;
       ierr = Fint(); CHKERRQ(ierr);
      for(int rr = 0;rr<row_mat;rr++) {
	if(RowGlob[rr].size()!=f_int[rr].size()) SETERRQ(PETSC_COMM_SELF,1,"data inconsistency");
	if(RowGlob[rr].size()==0) continue;
	ierr = VecSetValues(F_int,RowGlob[rr].size(),&(RowGlob[rr])[0],&(f_int[rr].data()[0]),ADD_VALUES); CHKERRQ(ierr);
      }
      
        PetscFunctionReturn(0);
    }


    virtual PetscErrorCode RhsAndLhs() {
      PetscFunctionBegin;
      //cout<<"RhsAndLhs "<<endl;
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
      //cout<<"Hi HI Hi"<<endl;
      
        
        //cout<<"F.size() "<<F.size()<<endl;
        //cout<<"F[0] "<<F[0]<<endl;
        
      ierr = OpStudentEnd(); CHKERRQ(ierr);
      PetscFunctionReturn(0);
    }

};

    
}

#endif //__THERMALFEMETHOD_HPP__
