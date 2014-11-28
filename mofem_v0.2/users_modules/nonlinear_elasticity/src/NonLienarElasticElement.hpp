/** \file ThermalElement.hpp 
 * \brief Operators and data structures for thermal analys
 *
 * Implementation of nonliear eleastic element.
 *
 */

/* Copyright (C) 2013, Lukasz Kaczmarczyk (likask AT wp.pl)
 * --------------------------------------------------------------
 *
 * This file is part of MoFEM.
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

#ifndef __NONLINEAR_ELASTIC_HPP
#define __NONLINEAR_ELASTIC_HPP

#ifndef WITH_ADOL_C
  #error "MoFEM need to be compiled with ADOL-C"
#endif 

/** \brief struture grouping operators and data used for calulation of mass (convective) element
  * \ingroup nonlinear_eleastic_elem
  *
  * In order to assemble matrices and right hand vectors, the loops over
  * elements, enetities over that elememnts and finally loop over intergration
  * points are executed.
  *
  * Following implementation separte those three cegories of loops and to eeach
  * loop attach operator.
  *
  */
struct NonlinearElasticElement {

  /// \brief  definition of volume element
  struct MyVolumeFE: public TetElementForcesAndSourcesCore {
    MyVolumeFE(FieldInterface &_mField): TetElementForcesAndSourcesCore(_mField) {}
    
    /** \brief it is used to calculate nb. of Gauss integartion points
     *
     * for more details pleas look 
     *   Reference:
     *
     * Albert Nijenhuis, Herbert Wilf,
     * Combinatorial Algorithms for Computers and Calculators,
     * Second Edition,
     * Academic Press, 1978,
     * ISBN: 0-12-519260-6,
     * LC: QA164.N54.
     *
     * More details about algorithm 
     * http://people.sc.fsu.edu/~jburkardt/cpp_src/gm_rule/gm_rule.html
    **/
    int getRule(int order) { return (order-1); };
  };
  
  MyVolumeFE feRhs; ///< cauclate right hand side for tetrahedral elements
  MyVolumeFE& getLoopFeRhs() { return feRhs; } ///< get rhs volume element 
  MyVolumeFE feLhs; //< calculate left hand side for tetrahedral elements
  MyVolumeFE& getLoopFeLhs() { return feLhs; } ///< get lhs volume element

  FieldInterface &mField;
  short int tAg;

  NonlinearElasticElement(
    FieldInterface &m_field,short int tag):
    feRhs(m_field),feLhs(m_field),
    mField(m_field),tAg(tag) {}

  /** \brief data for calulation het conductivity and heat capacity elements
    * \infroup mofem_forces_and_sources 
    */
  struct BlockData {
    double E;
    double PoissonRatio;
    Range tEts; ///< constatins elements in block set
  }; 
  map<int,BlockData> setOfBlocks; ///< maps block set id with appropiate BlockData

  /** \brief common data used by volume elements
    * \infroup mofem_forces_and_sources 
    */
  struct CommonData {
    map<string,vector<ublas::vector<double> > > dataAtGaussPts;
    map<string,vector<ublas::matrix<double> > > gradAtGaussPts;
    string spatialPositions;
    string meshPositions;
    vector<ublas::matrix<double> > P; ///< this need to be I Piola-Kirchoff stress tensor
    vector<vector<double*> > jacStressRowPtr;
    vector<ublas::matrix<double> > jacStress; ///< this is simply material tangent operator
  };
  CommonData commonData;

  struct OpGetDataAtGaussPts: public TetElementForcesAndSourcesCore::UserDataOperator {

    vector<ublas::vector<double> > &valuesAtGaussPts;
    vector<ublas::matrix<double> > &gradientAtGaussPts;
    const EntityType zeroAtType;

    OpGetDataAtGaussPts(const string field_name,
      vector<ublas::vector<double> > &values_at_gauss_pts,
      vector<ublas::matrix<double> > &gardient_at_gauss_pts):
      TetElementForcesAndSourcesCore::UserDataOperator(field_name),
      valuesAtGaussPts(values_at_gauss_pts),gradientAtGaussPts(gardient_at_gauss_pts),
      zeroAtType(MBVERTEX) {}

    /** \brief operator calulating deformation gradient
      *
      * temerature gradient is calculated multiplying direvatives of shape functions by degrees of freedom
      */
    PetscErrorCode doWork(
      int side,EntityType type,DataForcesAndSurcesCore::EntData &data) {
      PetscFunctionBegin;
      try {

 	int nb_dofs = data.getFieldData().size();
	if(nb_dofs == 0) {
	  PetscFunctionReturn(0);
	}
	int nb_gauss_pts = data.getN().size1();

	//initialize
	ublas::vector<double>& values = data.getFieldData();
	valuesAtGaussPts.resize(nb_gauss_pts);
	gradientAtGaussPts.resize(nb_gauss_pts);
	for(int gg = 0;gg<nb_gauss_pts;gg++) {
	  valuesAtGaussPts[gg].resize(3);
	  gradientAtGaussPts[gg].resize(3,3);
	}

        if(type == zeroAtType) { 
	  for(int gg = 0;gg<nb_gauss_pts;gg++) {
	    valuesAtGaussPts[gg].clear();
	    gradientAtGaussPts[gg].clear();
	  }
        }

	//cerr << valuesAtGaussPts[0] << " : ";
    
	for(int gg = 0;gg<nb_gauss_pts;gg++) {
	  ublas::vector<double> N = data.getN(gg,nb_dofs/3);
	  ublas::matrix<double> diffN = data.getDiffN(gg,nb_dofs/3);
	  for(int dd = 0;dd<nb_dofs/3;dd++) {
	    for(int rr1 = 0;rr1<3;rr1++) {
	      valuesAtGaussPts[gg][rr1] += N[dd]*values[3*dd+rr1];
	      for(int rr2 = 0;rr2<3;rr2++) {
		gradientAtGaussPts[gg](rr1,rr2) += diffN(dd,rr2)*values[3*dd+rr1];
	      }
	    }
	  }
	}

	//cerr << row_field_name << " " << col_field_name << endl;
	//cerr << side << " " << type << endl;
	//cerr << values << endl;
	//cerr << valuesAtGaussPts[0] << endl;

      } catch (const std::exception& ex) {
	ostringstream ss;
	ss << "throw in method: " << ex.what() << endl;
	SETERRQ(PETSC_COMM_SELF,1,ss.str().c_str());
      }

      PetscFunctionReturn(0);
    }

  };

  struct OpGetCommonDataAtGaussPts: public OpGetDataAtGaussPts {
    OpGetCommonDataAtGaussPts(const string field_name,CommonData &common_data):
      OpGetDataAtGaussPts(field_name,
      common_data.dataAtGaussPts[field_name], 
      common_data.gradAtGaussPts[field_name]) {}
  };
  
  struct FunctionsToCalulatePiolaKirchhoffI {

    double lambda,mu;
    ublas::matrix<adouble> F,C,E,S,P;

    PetscErrorCode CalulateC_CauchyDefromationTensor() {
      PetscFunctionBegin;
      C.resize(3,3);
      noalias(C) = prod(trans(F),F);
      PetscFunctionReturn(0);
    }

    PetscErrorCode CalulateE_GreenStrain() {
      PetscFunctionBegin;
      E.resize(3,3);
      noalias(E) = C;
      for(int dd = 0;dd<3;dd++) {
	E(dd,dd) -= 1;
      }
      E *= 0.5;
      PetscFunctionReturn(0);
    }

    //St. Venant–Kirchhoff Material
    PetscErrorCode CalculateS_PiolaKirchhoffII() {
      PetscFunctionBegin;
      adouble trE = 0;
      for(int dd = 0;dd<3;dd++) {
	trE += E(dd,dd);
      }
      S.resize(3,3);
      for(int dd = 0;dd<3;dd++) {
	S(dd,dd) = trE*lambda;
      }
      S += 2*mu*E;
      PetscFunctionReturn(0);
    }

    virtual PetscErrorCode CalualteP_PiolaKirchhoffI(
      const BlockData block_data,
      const NumeredMoFEMFiniteElement *fe_ptr) {
      PetscFunctionBegin;
      PetscErrorCode ierr;
      lambda = LAMBDA(block_data.E,block_data.PoissonRatio);
      mu = MU(block_data.E,block_data.PoissonRatio);
      //cerr << "lambda: " << lambda << endl;
      //cerr << "mu: " << mu << endl;
      ierr = CalulateC_CauchyDefromationTensor(); CHKERRQ(ierr);
      //cerr << "C: " << C << endl;
      ierr = CalulateE_GreenStrain(); CHKERRQ(ierr);
      //cerr << "E: " << E << endl;
      ierr = CalculateS_PiolaKirchhoffII(); CHKERRQ(ierr);
      //cerr << "S: " << S << endl;
      P.resize(3,3);
      noalias(P) = prod(F,S);
      //cerr << "P: " << P << endl;
      PetscFunctionReturn(0);
    }
  };

  struct OpJacobian: public TetElementForcesAndSourcesCore::UserDataOperator {

    BlockData &dAta;
    CommonData &commonData;
    FunctionsToCalulatePiolaKirchhoffI &fUn;
    int tAg;
    bool jAcobian;

    OpJacobian(
      const string field_name,
      BlockData &data,
      CommonData &common_data,
      FunctionsToCalulatePiolaKirchhoffI &fun,
      int tag,bool jacobian = true):
      TetElementForcesAndSourcesCore::UserDataOperator(field_name),
      dAta(data),commonData(common_data),fUn(fun),
      tAg(tag),jAcobian(jacobian) { }

    ublas::vector<double> active_varibles;

    PetscErrorCode doWork(
      int row_side,EntityType row_type,DataForcesAndSurcesCore::EntData &row_data) {
      PetscFunctionBegin;

      PetscErrorCode ierr;
      if(dAta.tEts.find(getMoFEMFEPtr()->get_ent()) == dAta.tEts.end()) {
	PetscFunctionReturn(0);
      }

      //do it only once, no need to repeat this for edges,faces or tets
      if(row_type != MBVERTEX) PetscFunctionReturn(0);

      int nb_dofs = row_data.getIndices().size();
      if(nb_dofs==0) PetscFunctionReturn(0);

      try {

	int nb_gauss_pts = row_data.getN().size1();
	commonData.P.resize(nb_gauss_pts);
	commonData.jacStressRowPtr.resize(nb_gauss_pts);
	commonData.jacStress.resize(nb_gauss_pts);
	int nb_active_variables = 0;

	for(int gg = 0;gg<nb_gauss_pts;gg++) {

	  if(gg == 0) {

	    //recorder on
	    trace_on(tAg);
	    
	    fUn.F.resize(3,3);
	    for(int dd1 = 0;dd1<3;dd1++) {
	      for(int dd2 = 0;dd2<3;dd2++) {
		fUn.F(dd1,dd2) <<= (commonData.gradAtGaussPts[commonData.spatialPositions][gg])(dd1,dd2);
		nb_active_variables++;
	      }
	    }
	    ierr = fUn.CalualteP_PiolaKirchhoffI(dAta,getMoFEMFEPtr()); CHKERRQ(ierr);
	    commonData.P[gg].resize(3,3);
	    for(int dd1 = 0;dd1<3;dd1++) {
	      for(int dd2 = 0;dd2<3;dd2++) {
		fUn.P(dd1,dd2) >>= (commonData.P[gg])(dd1,dd2);
	      }
	    }
	    
	    trace_off();
	    //recorder off

	  }

	  active_varibles.resize(nb_active_variables);
	  for(int dd1 = 0;dd1<3;dd1++) {
	    for(int dd2 = 0;dd2<3;dd2++) {
	      active_varibles(dd1*3+dd2) = (commonData.gradAtGaussPts[commonData.spatialPositions][gg])(dd1,dd2);
	    }
	  }

	  if(!jAcobian) {
	    if(gg>0) {
	      commonData.P[gg].resize(3,3);
	      int r;
	      //play recorder for values
	      r = function(tAg,9,nb_active_variables,&active_varibles[0],&commonData.P[gg](0,0));
	      if(r!=3) { // function is locally analytic
		SETERRQ1(PETSC_COMM_SELF,MOFEM_OPERATION_UNSUCCESSFUL,"ADOL-C function evaluation with error r = %d",r);
	      }
	    } 
	  } else {
	    commonData.jacStressRowPtr[gg].resize(9);
	    commonData.jacStress[gg].resize(9,nb_active_variables);
	    for(int dd = 0;dd<9;dd++) {
	      (commonData.jacStressRowPtr[gg])[dd] = &(commonData.jacStress[gg](dd,0));     
	    }
	    int r;
	    //play recorder for jacobians
	    r = jacobian(
	      tAg,9,nb_active_variables,
	      &active_varibles[0],&(commonData.jacStressRowPtr[gg])[0]);
	    if(r!=3) {
	      SETERRQ(PETSC_COMM_SELF,MOFEM_OPERATION_UNSUCCESSFUL,"ADOL-C function evaluation with error");
	    }
	  }

	}

      } catch (const std::exception& ex) {
	ostringstream ss;
	ss << "throw in method: " << ex.what() << endl;
	SETERRQ(PETSC_COMM_SELF,1,ss.str().c_str());
      }

      PetscFunctionReturn(0);
    }

  };

  struct OpRhs: public TetElementForcesAndSourcesCore::UserDataOperator {

    BlockData &dAta;
    CommonData &commonData;

    OpRhs(const string field_name,BlockData &data,CommonData &common_data):
      TetElementForcesAndSourcesCore::UserDataOperator(field_name),
      dAta(data),commonData(common_data) { }

    ublas::vector<double> nf;

    PetscErrorCode doWork(
      int row_side,EntityType row_type,DataForcesAndSurcesCore::EntData &row_data) {
      PetscFunctionBegin;

      PetscErrorCode ierr;
      if(dAta.tEts.find(getMoFEMFEPtr()->get_ent()) == dAta.tEts.end()) {
	PetscFunctionReturn(0);
      }
      if(row_data.getIndices().size()==0) PetscFunctionReturn(0);
      int nb_dofs = row_data.getIndices().size();

      try {

	nf.resize(nb_dofs);
	nf.clear();

	for(unsigned int gg = 0;gg<row_data.getN().size1();gg++) {
	  //diffN - on rows has degrees of freedom
	  //diffN - on columns has rerevatives direvatives of shape functin
	  const DataForcesAndSurcesCore::MatrixAdaptor &diffN = row_data.getDiffN(gg,nb_dofs/3);
	  const ublas::matrix<double>& P = commonData.P[gg];
	  //cerr << (commonData.gradAtGaussPts[commonData.spatialPositions][gg]) << endl;
	  //cerr << diffN << endl;
	  //cerr << P << endl;
	  double val = getVolume()*getGaussPts()(3,gg);
	  for(int dd = 0;dd<nb_dofs/3;dd++) {
	    for(int rr = 0;rr<3;rr++) {
	      for(int nn = 0;nn<3;nn++) {
		nf[3*dd+rr] += val*diffN(dd,nn)*P(rr,nn);
	      }
	    }
	  }
	}

	if((unsigned int)nb_dofs > 3*row_data.getN().size2()) {
	  SETERRQ(PETSC_COMM_SELF,1,"data inconsistency");
	}
	ierr = VecSetValues(getFEMethod()->snes_f,nb_dofs,
	  &row_data.getIndices()[0],&nf[0],ADD_VALUES); CHKERRQ(ierr);

      } catch (const std::exception& ex) {
	ostringstream ss;
	ss << "throw in method: " << ex.what() << endl;
	SETERRQ(PETSC_COMM_SELF,1,ss.str().c_str());
      }

      PetscFunctionReturn(0);
    }

  };

  struct OpLhs_dx: public TetElementForcesAndSourcesCore::UserDataOperator {

    BlockData &dAta;
    CommonData &commonData;
    int tAg;

    OpLhs_dx(const string vel_field,const string field_name,BlockData &data,CommonData &common_data):
      TetElementForcesAndSourcesCore::UserDataOperator(vel_field,field_name),
      dAta(data),commonData(common_data) { /*symm = false;*/  }

    ublas::matrix<double> k,trans_k,jac,F;

    virtual PetscErrorCode getJac(DataForcesAndSurcesCore::EntData &col_data,int gg) {
      PetscFunctionBegin;
      int nb_col = col_data.getFieldData().size();
      const DataForcesAndSurcesCore::MatrixAdaptor diffN = col_data.getDiffN(gg,nb_col/3);
      F.resize(3,3);
      for(int dd = 0;dd<nb_col/3;dd++) {
	  F.clear();
	  for(int rr = 0;rr<3;rr++) {
	    F(rr,0) = diffN(dd,0);
	    F(rr,1) = diffN(dd,1);
	    F(rr,2) = diffN(dd,2);
	    for(int ii = 0;ii<9;ii++) {
	      for(int jj = 0;jj<9;jj++) {
		//This project dirvative \frac{\partial P}{\partial F}, that is
		//\frac{\partial P}{\partial x_DOF} =  \frac{\partial P}{\partial F}\frac{\partial F}{\partial x_DOF},
		//where second therm \frac{\partial F}{\partial x_DOF} is derivative of shape function
		jac(ii,3*dd+rr) += commonData.jacStress[gg](ii,jj)*F.data()[jj];
	      }
	    }	
	  }
	}
      PetscFunctionReturn(0);
    }

    PetscErrorCode doWork(
      int row_side,int col_side,
      EntityType row_type,EntityType col_type,
      DataForcesAndSurcesCore::EntData &row_data,
      DataForcesAndSurcesCore::EntData &col_data) {      
      PetscFunctionBegin;

      PetscErrorCode ierr;

      if(dAta.tEts.find(getMoFEMFEPtr()->get_ent()) == dAta.tEts.end()) {
	PetscFunctionReturn(0);
      }
      
      int nb_row = row_data.getIndices().size();
      int nb_col = col_data.getIndices().size();
      if(nb_row == 0) PetscFunctionReturn(0);
      if(nb_col == 0) PetscFunctionReturn(0);

      try {

	k.resize(nb_row,nb_col);
	k.clear();
	jac.resize(9,nb_col);

	for(unsigned int gg = 0;gg<row_data.getN().size1();gg++) {

	  ierr = getJac(col_data,gg); CHKERRQ(ierr);
	  double val = getVolume()*getGaussPts()(3,gg);
	  jac *= val;

	  const DataForcesAndSurcesCore::MatrixAdaptor &diffN = row_data.getDiffN(gg,nb_row/3);

	  { //integrate element stiffnes matrix
	    for(int dd1 = 0;dd1<nb_row/3;dd1++) {
	      for(int rr1 = 0;rr1<3;rr1++) {
		for(int dd2 = 0;dd2<nb_col/3;dd2++) {
		  for(int rr2 = 0;rr2<3;rr2++) {
		    k(3*dd1+rr1,3*dd2+rr2) += 
		      diffN(dd1,0)*jac(3*rr1+0,3*dd2+rr2)+
		      diffN(dd1,1)*jac(3*rr1+1,3*dd2+rr2)+
   		      diffN(dd1,2)*jac(3*rr1+2,3*dd2+rr2);
		  }
		}
	      }
	    }
	  }

	}

	ierr = MatSetValues(getFEMethod()->snes_B,
	  nb_row,&row_data.getIndices()[0],
	  nb_col,&col_data.getIndices()[0],
	  &k(0,0),ADD_VALUES); CHKERRQ(ierr);
	//is symmetric
        if(row_side != col_side || row_type != col_type) {
          trans_k.resize(nb_col,nb_row);
          noalias(trans_k) = trans( k );
          ierr = MatSetValues(
                 getFEMethod()->snes_B,
                 nb_col,&col_data.getIndices()[0],
                 nb_row,&row_data.getIndices()[0],
                 &trans_k(0,0),ADD_VALUES); CHKERRQ(ierr);
        }


      } catch (const std::exception& ex) {
	ostringstream ss;
	ss << "throw in method: " << ex.what() << endl;
	SETERRQ(PETSC_COMM_SELF,1,ss.str().c_str());
      }

      PetscFunctionReturn(0);
    }

  };

  PetscErrorCode setBlocks() {
    PetscFunctionBegin;
    ErrorCode rval;
    PetscErrorCode ierr;
  
    for(_IT_CUBITMESHSETS_BY_BCDATA_TYPE_FOR_LOOP_(mField,BLOCKSET|MAT_ELASTICSET,it)) {
      Mat_Elastic mydata;
      ierr = it->get_attribute_data_structure(mydata); CHKERRQ(ierr);
      int id = it->get_msId();
      EntityHandle meshset = it->get_meshset();
      rval = mField.get_moab().get_entities_by_type(meshset,MBTET,setOfBlocks[id].tEts,true); CHKERR_PETSC(rval);
      setOfBlocks[id].E = mydata.data.Young;
      setOfBlocks[id].PoissonRatio = mydata.data.Poisson;
      //cerr << setOfBlocks[id].tEts << endl;
    }


    PetscFunctionReturn(0);
  }

  PetscErrorCode addElement(string element_name,
    string spatial_position_field_name,
    string material_position_field_name = "MESH_NODE_POSITIONS",bool ale = false) {
    PetscFunctionBegin;

    PetscErrorCode ierr;
    //ErrorCode rval;

    ierr = mField.add_finite_element(element_name,MF_ZERO); CHKERRQ(ierr);
    ierr = mField.modify_finite_element_add_field_row(element_name,spatial_position_field_name); CHKERRQ(ierr);
    ierr = mField.modify_finite_element_add_field_col(element_name,spatial_position_field_name); CHKERRQ(ierr);
    ierr = mField.modify_finite_element_add_field_data(element_name,spatial_position_field_name); CHKERRQ(ierr);
    if(mField.check_field(material_position_field_name)) {
      if(ale) {
	ierr = mField.modify_finite_element_add_field_row(element_name,material_position_field_name); CHKERRQ(ierr);
	ierr = mField.modify_finite_element_add_field_col(element_name,material_position_field_name); CHKERRQ(ierr);
      }
      ierr = mField.modify_finite_element_add_field_data(element_name,material_position_field_name); CHKERRQ(ierr);
    }

    map<int,BlockData>::iterator sit = setOfBlocks.begin();
    for(;sit!=setOfBlocks.end();sit++) {
      ierr = mField.add_ents_to_finite_element_by_TETs(sit->second.tEts,element_name); CHKERRQ(ierr);
    }

    PetscFunctionReturn(0);
  }

  PetscErrorCode setOperators(
    FunctionsToCalulatePiolaKirchhoffI &fun,
    string spatial_position_field_name,
    string material_position_field_name = "MESH_NODE_POSITIONS",
    bool ale = false) {
    PetscFunctionBegin;

    commonData.spatialPositions = spatial_position_field_name;
    commonData.meshPositions = material_position_field_name;

    //Rhs
    feRhs.get_op_to_do_Rhs().push_back(new OpGetCommonDataAtGaussPts(spatial_position_field_name,commonData));
    if(mField.check_field(material_position_field_name)) {
      feRhs.get_op_to_do_Rhs().push_back(new OpGetCommonDataAtGaussPts(material_position_field_name,commonData));
    }
    map<int,BlockData>::iterator sit = setOfBlocks.begin();
    for(;sit!=setOfBlocks.end();sit++) {
      feRhs.get_op_to_do_Rhs().push_back(new OpJacobian(spatial_position_field_name,sit->second,commonData,fun,tAg,false));
      feRhs.get_op_to_do_Rhs().push_back(new OpRhs(spatial_position_field_name,sit->second,commonData));
    }

    //Lhs
    feLhs.get_op_to_do_Rhs().push_back(new OpGetCommonDataAtGaussPts(spatial_position_field_name,commonData));
    if(mField.check_field(material_position_field_name)) {
      feLhs.get_op_to_do_Rhs().push_back(new OpGetCommonDataAtGaussPts(material_position_field_name,commonData));
    }
    sit = setOfBlocks.begin();
    for(;sit!=setOfBlocks.end();sit++) {
      feLhs.get_op_to_do_Rhs().push_back(new OpJacobian(spatial_position_field_name,sit->second,commonData,fun,tAg));
      feLhs.get_op_to_do_Lhs().push_back(new OpLhs_dx(spatial_position_field_name,spatial_position_field_name,sit->second,commonData));
    }

    PetscFunctionReturn(0);
  }

};


#endif //__NONLINEAR_ELASTIC_HPP

/***************************************************************************//**
 * \defgroup nonlinear_eleastic_elem Non-Linear Elastic Element 
 * \ingroup mofem_forces_and_sources 
 ******************************************************************************/



