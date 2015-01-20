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

namespace ObosleteUsersModules {

struct ElasticFEMethod: public FEMethod_UpLevelStudent {

    FieldInterface& mField;
    bool propeties_from_BLOCKSET_MAT_ELASTICSET;
    double lambda,mu;
    string fieldName;

    ElasticFEMethod(FieldInterface& _mField,Mat _Aij,Vec _X,Vec _F,double _lambda,double _mu, string _field_name = "DISPLACEMENT"):
      FEMethod_UpLevelStudent(_mField.get_moab(),1),mField(_mField),lambda(_lambda),mu(_mu),fieldName(_field_name) {

      snes_B = _Aij;
      snes_x = _X;
      snes_f = _F;

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

      if(snes_f!=PETSC_NULL) {
	//VEC & MAT Options
	//If index is set to -1 ingonre its assembly
	VecSetOption(snes_f, VEC_IGNORE_NEGATIVE_INDICES,PETSC_TRUE); 
      }

      propeties_from_BLOCKSET_MAT_ELASTICSET = false;
      for(_IT_CUBITMESHSETS_BY_BCDATA_TYPE_FOR_LOOP_(mField,BLOCKSET|MAT_ELASTICSET,it)) {
	propeties_from_BLOCKSET_MAT_ELASTICSET = true;
      }

    }; 

    ErrorCode rval;
    
    ParallelComm* pcomm;
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

    vector<DofIdx> DirichletBC;

    vector<double> g_NTET;
    const double *G_TET_W;

    virtual PetscErrorCode calculateD(double _lambda,double _mu) {
      PetscFunctionBegin;

      D = _lambda*D_lambda + _mu*D_mu;
      //cerr << D_lambda << endl;
      //cerr << D_mu << endl;
      //cerr << D << endl;

      PetscFunctionReturn(0);
    }

    virtual PetscErrorCode GetMatParameters(double *_lambda,double *_mu,
      double *_User1 = NULL,double *_User2 = NULL) {
      PetscFunctionBegin;

      *_lambda = lambda;
      *_mu = mu;


      if(propeties_from_BLOCKSET_MAT_ELASTICSET) {
	EntityHandle ent = fePtr->get_ent();
	for(_IT_CUBITMESHSETS_BY_BCDATA_TYPE_FOR_LOOP_(mField,BLOCKSET|MAT_ELASTICSET,it)) {

	  Mat_Elastic mydata;
	  ierr = it->get_attribute_data_structure(mydata); CHKERRQ(ierr);

	  Range meshsets;
	  rval = moab.get_entities_by_type(it->meshset,MBENTITYSET,meshsets,true); CHKERR_PETSC(rval);
	  meshsets.insert(it->meshset);
	  for(Range::iterator mit = meshsets.begin();mit != meshsets.end(); mit++) {
	    if( moab.contains_entities(*mit,&ent,1) ) {
	      *_lambda = LAMBDA(mydata.data.Young,mydata.data.Poisson);
	      *_mu = MU(mydata.data.Young,mydata.data.Poisson);
	      if(_User1!=NULL) *_User1 = mydata.data.User1;
	      if(_User2!=NULL) *_User2 = mydata.data.User2;
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

    virtual PetscErrorCode GetMatricesRows() {
      PetscFunctionBegin;
      //indicies ROWS
      row_mat = 0;
      ierr = GetRowGlobalIndices(fieldName,RowGlob[row_mat]); CHKERRQ(ierr);
      ierr = GetRowLocalIndices(fieldName,RowLocal[row_mat]); CHKERRQ(ierr);
      ierr = GetGaussRowNMatrix(fieldName,rowNMatrices[row_mat]); CHKERRQ(ierr);
      ierr = GetGaussRowDiffNMatrix(fieldName,rowDiffNMatrices[row_mat]); CHKERRQ(ierr);
      //HO gemometry
      ierr = GetHierarchicalGeometryApproximation_ApplyToDiffShapeFunction(3,invH,rowDiffNMatrices[row_mat]); CHKERRQ(ierr);
      //
      ierr = MakeBMatrix3D(fieldName,rowDiffNMatrices[row_mat],rowBMatrices[row_mat]);  CHKERRQ(ierr);
      row_mat++;
      for(int ee = 0;ee<6;ee++) { //edges matrices
	RowGlob[row_mat].resize(0);
	ierr = GetRowGlobalIndices(fieldName,MBEDGE,RowGlob[row_mat],ee); CHKERRQ(ierr);
        RowLocal[row_mat].resize(0);
	ierr = GetRowLocalIndices(fieldName,MBEDGE,RowLocal[row_mat],ee); CHKERRQ(ierr);
	if(RowGlob[row_mat].size()!=0) {
	  ierr = GetGaussRowNMatrix(fieldName,MBEDGE,rowNMatrices[row_mat],ee); CHKERRQ(ierr);
	  ierr = GetGaussRowDiffNMatrix(fieldName,MBEDGE,rowDiffNMatrices[row_mat],ee); CHKERRQ(ierr);
	  //HO gemometry
	  ierr = GetHierarchicalGeometryApproximation_ApplyToDiffShapeFunction(3,invH,rowDiffNMatrices[row_mat]); CHKERRQ(ierr);
	  //
	  ierr = MakeBMatrix3D(fieldName,rowDiffNMatrices[row_mat],rowBMatrices[row_mat]);  CHKERRQ(ierr);
	  //cerr << rowDiffNMatrices[row_mat][0] << endl;
	}
	row_mat++;
      }
      for(int ff = 0;ff<4;ff++) { //faces matrices
	RowGlob[row_mat].resize(0);
	ierr = GetRowGlobalIndices(fieldName,MBTRI,RowGlob[row_mat],ff); CHKERRQ(ierr);
	RowLocal[row_mat].resize(0);
	ierr = GetRowLocalIndices(fieldName,MBTRI,RowLocal[row_mat],ff); CHKERRQ(ierr);
	if(RowGlob[row_mat].size()!=0) {
	  ierr = GetGaussRowNMatrix(fieldName,MBTRI,rowNMatrices[row_mat],ff); CHKERRQ(ierr);
	  ierr = GetGaussRowDiffNMatrix(fieldName,MBTRI,rowDiffNMatrices[row_mat],ff); CHKERRQ(ierr);
	  //HO gemometry
	  ierr = GetHierarchicalGeometryApproximation_ApplyToDiffShapeFunction(3,invH,rowDiffNMatrices[row_mat]); CHKERRQ(ierr);
	  //
	  ierr = MakeBMatrix3D(fieldName,rowDiffNMatrices[row_mat],rowBMatrices[row_mat]);  CHKERRQ(ierr);
	}
	row_mat++;
      }
      RowGlob[row_mat].resize(0);
      ierr = GetRowGlobalIndices(fieldName,MBTET,RowGlob[row_mat]); CHKERRQ(ierr);
      RowLocal[row_mat].resize(0);
      ierr = GetRowLocalIndices(fieldName,MBTET,RowLocal[row_mat]); CHKERRQ(ierr);
      if(RowGlob[row_mat].size() != 0) { //volume matrices
	ierr = GetGaussRowNMatrix(fieldName,MBTET,rowNMatrices[row_mat]); CHKERRQ(ierr);
	ierr = GetGaussRowDiffNMatrix(fieldName,MBTET,rowDiffNMatrices[row_mat]); CHKERRQ(ierr);
	//HO gemometry
	ierr = GetHierarchicalGeometryApproximation_ApplyToDiffShapeFunction(3,invH,rowDiffNMatrices[row_mat]); CHKERRQ(ierr);
	//
	ierr = MakeBMatrix3D(fieldName,rowDiffNMatrices[row_mat],rowBMatrices[row_mat]);  CHKERRQ(ierr);
      }
      row_mat++;
      PetscFunctionReturn(0);
    }

    virtual PetscErrorCode GetMatricesCols() {
      PetscFunctionBegin;
      //Higher order approximation of geometry
      //ierr = GetHierarchicalGeometryApproximation(invH,detH); CHKERRQ(ierr);
      //indicies COLS
      col_mat = 0;
      ierr = GetColGlobalIndices(fieldName,ColGlob[col_mat]); CHKERRQ(ierr);
      ierr = GetColLocalIndices(fieldName,ColLocal[col_mat]); CHKERRQ(ierr);
      ierr = GetGaussColNMatrix(fieldName,colNMatrices[col_mat]); CHKERRQ(ierr);
      ierr = GetGaussColDiffNMatrix(fieldName,colDiffNMatrices[col_mat]); CHKERRQ(ierr);
      //HO gemometry
      ierr = GetHierarchicalGeometryApproximation_ApplyToDiffShapeFunction(3,invH,colDiffNMatrices[col_mat]); CHKERRQ(ierr);
      //
      ierr = MakeBMatrix3D(fieldName,colDiffNMatrices[col_mat],colBMatrices[col_mat]);  CHKERRQ(ierr);
      col_mat++;
      for(int ee = 0;ee<6;ee++) { //edges matrices
	ierr = GetColGlobalIndices(fieldName,MBEDGE,ColGlob[col_mat],ee); CHKERRQ(ierr);
	ierr = GetColLocalIndices(fieldName,MBEDGE,ColLocal[col_mat],ee); CHKERRQ(ierr);
	if(ColGlob[col_mat].size()!=0) {
	  ierr = GetGaussColNMatrix(fieldName,MBEDGE,colNMatrices[col_mat],ee); CHKERRQ(ierr);
	  ierr = GetGaussColDiffNMatrix(fieldName,MBEDGE,colDiffNMatrices[col_mat],ee); CHKERRQ(ierr);
	  //HO gemometry
	  ierr = GetHierarchicalGeometryApproximation_ApplyToDiffShapeFunction(3,invH,colDiffNMatrices[col_mat]); CHKERRQ(ierr);
	  //
	  ierr = MakeBMatrix3D(fieldName,colDiffNMatrices[col_mat],colBMatrices[col_mat]);  CHKERRQ(ierr);
	}
	col_mat++;
      }
      for(int ff = 0;ff<4;ff++) { //faces matrices
	ierr = GetColGlobalIndices(fieldName,MBTRI,ColGlob[col_mat],ff); CHKERRQ(ierr);
	ierr = GetColLocalIndices(fieldName,MBTRI,ColLocal[col_mat],ff); CHKERRQ(ierr);
	if(ColGlob[col_mat].size()!=0) {
	  ierr = GetGaussColNMatrix(fieldName,MBTRI,colNMatrices[col_mat],ff); CHKERRQ(ierr);
	  ierr = GetGaussColDiffNMatrix(fieldName,MBTRI,colDiffNMatrices[col_mat],ff); CHKERRQ(ierr);
	  //HO gemometry
	  ierr = GetHierarchicalGeometryApproximation_ApplyToDiffShapeFunction(3,invH,colDiffNMatrices[col_mat]); CHKERRQ(ierr);
	  //
	  ierr = MakeBMatrix3D(fieldName,colDiffNMatrices[col_mat],colBMatrices[col_mat]);  CHKERRQ(ierr);
	}
	col_mat++;
      }
      ierr = GetColGlobalIndices(fieldName,MBTET,ColGlob[col_mat]); CHKERRQ(ierr);
      ierr = GetColLocalIndices(fieldName,MBTET,ColLocal[col_mat]); CHKERRQ(ierr);
      if(ColGlob[col_mat].size() != 0) { //volume matrices
	ierr = GetGaussColNMatrix(fieldName,MBTET,colNMatrices[col_mat]); CHKERRQ(ierr);
	ierr = GetGaussColDiffNMatrix(fieldName,MBTET,colDiffNMatrices[col_mat]); CHKERRQ(ierr);
	//HO gemometry
	ierr = GetHierarchicalGeometryApproximation_ApplyToDiffShapeFunction(3,invH,colDiffNMatrices[col_mat]); CHKERRQ(ierr);
	//
	ierr = MakeBMatrix3D(fieldName,colDiffNMatrices[col_mat],colBMatrices[col_mat]);  CHKERRQ(ierr);
      }
      col_mat++;

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
      ierr = calculateD(_lambda,_mu); CHKERRQ(ierr);

      K.resize(row_mat,col_mat);
      int g_dim = g_NTET.size()/4;
      for(int rr = 0;rr<row_mat;rr++) {
	if(RowGlob[rr].size()==0) continue;
	for(int gg = 0;gg<g_dim;gg++) {
	  ublas::matrix<FieldData> &row_Mat = (rowBMatrices[rr])[gg];
	  double w = V*G_TET_W[gg];
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
	    if(ColGlob[cc].size()==0) continue;
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
	if(RowGlob[rr].size()==0) continue;
	for(int cc = rr;cc<col_mat;cc++) {
	  if(ColGlob[cc].size()==0) continue;
	  if(RowGlob[rr].size()!=K(rr,cc).size1()) SETERRQ(PETSC_COMM_SELF,1,"data inconsistency");
	  if(ColGlob[cc].size()!=K(rr,cc).size2()) SETERRQ(PETSC_COMM_SELF,1,"data inconsistency");
	  ierr = MatSetValues(snes_B,RowGlob[rr].size(),&(RowGlob[rr])[0],ColGlob[cc].size(),&(ColGlob[cc])[0],&(K(rr,cc).data())[0],ADD_VALUES); CHKERRQ(ierr);
	  if(rr!=cc) {
	    K(cc,rr) = trans(K(rr,cc));
	    ierr = MatSetValues(snes_B,ColGlob[cc].size(),&(ColGlob[cc])[0],RowGlob[rr].size(),&(RowGlob[rr])[0],&(K(cc,rr).data())[0],ADD_VALUES); CHKERRQ(ierr);
	  }
	}
      }
      PetscFunctionReturn(0);
    }

    virtual PetscErrorCode Rhs() {
      PetscFunctionBegin;
      ierr = Fint(snes_f); CHKERRQ(ierr);
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
      ierr = calculateD(_lambda,_mu); CHKERRQ(ierr);

      //Gradient at Gauss points; 
      vector< ublas::matrix< FieldData > > GradU_at_GaussPt;
      ierr = GetGaussDiffDataVector(fieldName,GradU_at_GaussPt); CHKERRQ(ierr);
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
	  double w = V*G_TET_W[gg];
	  ublas::vector<FieldData> VoightStress = prod(w*D,VoightStrain);
	  //BT * VoigtStress
	  f_int.resize(row_mat);
	  for(int rr = 0;rr<row_mat;rr++) {
	    if(RowGlob[rr].size()==0) continue;
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
	  if(RowGlob[rr].size()==0) continue;
	  if(RowGlob[rr].size()!=f_int[rr].size()) SETERRQ(PETSC_COMM_SELF,1,"data inconsistency");
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

    ublas::matrix<double> gaussPts;
    virtual PetscErrorCode Get_g_NTET() {
      PetscFunctionBegin;

      int order = 1;
      for(_IT_GET_FEDATA_BY_NAME_DOFS_FOR_LOOP_(this,fieldName,dof)) {
	order = max(order,dof->get_max_order());
      }

      int rule = max(0,order-1);
      if( 2*rule + 1 < 2*(order-1) ) {
	SETERRQ2(PETSC_COMM_SELF,1,"wrong rule %d %d",order,rule);
      }
      int nb_gauss_pts = gm_rule_size(rule,3);
      if(gaussPts.size2() == (unsigned int)nb_gauss_pts) {
	PetscFunctionReturn(0);
      }
      gaussPts.resize(4,nb_gauss_pts);
      ierr = Grundmann_Moeller_integration_points_3D_TET(
	rule,&gaussPts(0,0),&gaussPts(1,0),&gaussPts(2,0),&gaussPts(3,0)); CHKERRQ(ierr);

      g_NTET.resize(4*nb_gauss_pts);
      ierr = ShapeMBTET(&g_NTET[0],&gaussPts(0,0),&gaussPts(1,0),&gaussPts(2,0),nb_gauss_pts); CHKERRQ(ierr);
      G_TET_W = &gaussPts(3,0);

      PetscFunctionReturn(0);
    }

    PetscErrorCode preProcess() {
      PetscFunctionBegin;

      g_NTET.resize(4*45);
      ierr = ShapeMBTET(&g_NTET[0],G_TET_X45,G_TET_Y45,G_TET_Z45,45); CHKERRQ(ierr);
      G_TET_W = G_TET_W45;
			
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

      PetscFunctionReturn(0);
    }

    PetscErrorCode postProcess() {
      PetscFunctionBegin;

      switch(snes_ctx) {
        case CTX_SNESNONE: {
	  // Note MAT_FLUSH_ASSEMBLY
	  ierr = MatAssemblyBegin(snes_B,MAT_FLUSH_ASSEMBLY); CHKERRQ(ierr);
	  ierr = MatAssemblyEnd(snes_B,MAT_FLUSH_ASSEMBLY); CHKERRQ(ierr);
	  ierr = VecAssemblyBegin(snes_f); CHKERRQ(ierr);
	  ierr = VecAssemblyEnd(snes_f); CHKERRQ(ierr);
        }
	break;
        case CTX_SNESSETFUNCTION: {
	  ierr = VecAssemblyBegin(snes_f); CHKERRQ(ierr);
	  ierr = VecAssemblyEnd(snes_f); CHKERRQ(ierr);
	}
	break;
        case CTX_SNESSETJACOBIAN: {
	  // Note MAT_FLUSH_ASSEMBLY
	  ierr = MatAssemblyBegin(snes_B,MAT_FLUSH_ASSEMBLY); CHKERRQ(ierr);
	  ierr = MatAssemblyEnd(snes_B,MAT_FLUSH_ASSEMBLY); CHKERRQ(ierr);
	}
        break;
        default:
          SETERRQ(PETSC_COMM_SELF,1,"not implemented");
      }

      PetscFunctionReturn(0);
    }

    PetscErrorCode operator()() {
      PetscFunctionBegin;

      ierr = Get_g_NTET(); CHKERRQ(ierr);
      ierr = OpStudentStart_TET(g_NTET); CHKERRQ(ierr);
      ierr = GetMatrices(); CHKERRQ(ierr);

      switch(snes_ctx) {
        case CTX_SNESNONE: {
	  ierr = RhsAndLhs(); CHKERRQ(ierr);
        }
        case CTX_SNESSETFUNCTION: {
	  ierr = Rhs(); CHKERRQ(ierr);
	}
	break;
        case CTX_SNESSETJACOBIAN: {
	  ierr = Lhs(); CHKERRQ(ierr);
	}
        break;
        default:
          SETERRQ(PETSC_COMM_SELF,1,"not implemented");
      }

      ierr = OpStudentEnd(); CHKERRQ(ierr);
      PetscFunctionReturn(0);
    }

};

    
}

#endif //__ELASTICFEMETHOD_HPP__
