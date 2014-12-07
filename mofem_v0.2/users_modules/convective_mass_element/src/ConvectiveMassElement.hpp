/** \file ThermalElement.hpp 
 * \brief Operators and data structures for thermal analys
 *
 * Implementation of convective mass element
 *
 */

/* Copyright (C) 2014, Lukasz Kaczmarczyk (likask AT wp.pl)
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

#ifndef __CONVECTIVE_MASS_ELEMENT_HPP
#define __CONVECTIVE_MASS_ELEMENT_HPP

#ifndef WITH_ADOL_C
  #error "MoFEM need to be compiled with ADOL-C"
#endif 

/** \brief struture grouping operators and data used for calulation of mass (convective) element
  * \ingroup convective_mass_elem
  *
  * In order to assemble matrices and right hand vectors, the loops over
  * elements, enetities over that elememnts and finally loop over intergration
  * points are executed.
  *
  * Following implementation separte those three cegories of loops and to eeach
  * loop attach operator.
  *
  */
struct ConvectiveMassElement {

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
    int getRule(int order) { return order; };
  };
  
  MyVolumeFE feMassRhs; ///< cauclate right hand side for tetrahedral elements
  MyVolumeFE& getLoopFeMassRhs() { return feMassRhs; } ///< get rhs volume element 
  MyVolumeFE feMassLhs; //< calculate left hand side for tetrahedral elements
  MyVolumeFE& getLoopFeMassLhs() { return feMassLhs; } ///< get lhs volume element

  MyVolumeFE feVelRhs; ///< cauclate right hand side for tetrahedral elements
  MyVolumeFE& getLoopFeVelRhs() { return feVelRhs; } ///< get rhs volume element 
  MyVolumeFE feVelLhs; //< calculate left hand side for tetrahedral elements
  MyVolumeFE& getLoopFeVelLhs() { return feVelLhs; } ///< get lhs volume element

  MyVolumeFE feTRhs; ///< cauclate right hand side for tetrahedral elements
  MyVolumeFE& getLoopFeTRhs() { return feTRhs; } ///< get rhs volume element 
  MyVolumeFE feTLhs; //< calculate left hand side for tetrahedral elements
  MyVolumeFE& getLoopFeTLhs() { return feTLhs; } ///< get lhs volume element

  FieldInterface &mField;
  short int tAg;

  ConvectiveMassElement(
    FieldInterface &m_field,short int tag):
    feMassRhs(m_field),feMassLhs(m_field),
    feVelRhs(m_field),feVelLhs(m_field),
    feTRhs(m_field),feTLhs(m_field),
    mField(m_field),tAg(tag) {}

  /** \brief data for calulation het conductivity and heat capacity elements
    * \infroup mofem_forces_and_sources 
    */
  struct BlockData {
    double rho0; ///< reference density
    ublas::vector<double> a0; //< constant acceleration
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
    string spatialVelocities;
    vector<ublas::vector<double> > valVel;
    vector<vector<double*> > jacVelRowPtr;
    vector<ublas::matrix<double> > jacVel;
    vector<ublas::vector<double> > valMass;
    vector<vector<double*> > jacMassRowPtr;
    vector<ublas::matrix<double> > jacMass;
    vector<double > valT;
    vector<ublas::vector<double> > jacEshelby;
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
  
  struct CommonFunctions {

    template<typename TYPE> 
    PetscErrorCode dEterminatnt(ublas::matrix<TYPE> a,TYPE &det) {
      PetscFunctionBegin;
      //a11a22a33
      //+a21a32a13
      //+a31a12a23
      //-a11a32a23
      //-a31a22a13
      //-a21a12a33
      //http://www.cg.info.hiroshima-cu.ac.jp/~miyazaki/knowledge/teche23.html
      //http://mathworld.wolfram.com/MatrixInverse.html
      det = a(0,0)*a(1,1)*a(2,2)
        +a(1,0)*a(2,1)*a(0,2)
        +a(2,0)*a(0,1)*a(1,2)
        -a(0,0)*a(2,1)*a(1,2)
        -a(2,0)*a(1,1)*a(0,2)
        -a(1,0)*a(0,1)*a(2,2);
      PetscFunctionReturn(0);
    }
  
    template<typename TYPE> 
    PetscErrorCode iNvert(TYPE det,ublas::matrix<TYPE> a,ublas::matrix<TYPE> &inv_a) {
      PetscFunctionBegin;
      //PetscErrorCode ierr;
      inv_a.resize(3,3);
      //http://www.cg.info.hiroshima-cu.ac.jp/~miyazaki/knowledge/teche23.html
      //http://mathworld.wolfram.com/MatrixInverse.html
      inv_a(0,0) = a(1,1)*a(2,2)-a(1,2)*a(2,1);
      inv_a(0,1) = a(0,2)*a(2,1)-a(0,1)*a(2,2);
      inv_a(0,2) = a(0,1)*a(1,2)-a(0,2)*a(1,1);
      inv_a(1,0) = a(1,2)*a(2,0)-a(1,0)*a(2,2);
      inv_a(1,1) = a(0,0)*a(2,2)-a(0,2)*a(2,0);
      inv_a(1,2) = a(0,2)*a(1,0)-a(0,0)*a(1,2);
      inv_a(2,0) = a(1,0)*a(2,1)-a(1,1)*a(2,0);
      inv_a(2,1) = a(0,1)*a(2,0)-a(0,0)*a(2,1);
      inv_a(2,2) = a(0,0)*a(1,1)-a(0,1)*a(1,0);
      //TYPE det;
      //ierr = dEterminatnt(a,det); CHKERRQ(ierr);
      inv_a /= det;
      PetscFunctionReturn(0);
    }
  
    template<typename TYPE> 
    PetscErrorCode calculateMomentumRate(
      double rho0,ublas::vector<double>& a0,
      ublas::vector<TYPE>& a,ublas::matrix<TYPE>& grad_v,
      ublas::vector<TYPE>& dot_W,
      ublas::matrix<TYPE>& H,ublas::matrix<TYPE>& invH,
      ublas::matrix<TYPE>& h,
      ublas::matrix<TYPE>& F,
      ublas::vector<TYPE>& dp_dt) {
      PetscFunctionBegin;
      PetscErrorCode ierr;
      //calulate gradient of deformation
      noalias(F) = prod(h,invH);
      TYPE detF;
      ierr = dEterminatnt(F,detF); CHKERRQ(ierr);
      //calulate current density
      TYPE rho = rho0*detF*rho0;
      //momentum rate
      noalias(dp_dt) = rho*(a0 + a - prod(grad_v,dot_W));
      PetscFunctionReturn(0);
    }
  
    template<typename TYPE> 
    PetscErrorCode calculateFuncUnderIntegral(
      double rho0,ublas::vector<double>& a0,
      ublas::vector<TYPE>& a,
      ublas::matrix<TYPE>& grad_v,
      ublas::vector<TYPE>& dot_W,
      ublas::matrix<TYPE>& H,
      TYPE& detH,
      ublas::matrix<TYPE>& invH,
      ublas::matrix<TYPE>& h,
      ublas::matrix<TYPE>& F,
      ublas::vector<TYPE>& dp_dt,
      ublas::vector<TYPE>& f) {
      PetscFunctionBegin;
      PetscErrorCode ierr;
      ierr = calculateMomentumRate(
        rho0,a0,a,grad_v,dot_W,H,invH,h,F,dp_dt); CHKERRQ(ierr);
      noalias(f) = dp_dt*detH;
      PetscFunctionReturn(0);
    }

  };

  struct OpMassJacobian: public TetElementForcesAndSourcesCore::UserDataOperator,CommonFunctions {

    BlockData &dAta;
    CommonData &commonData;
    int tAg;
    bool jAcobian;

    OpMassJacobian(const string field_name,BlockData &data,CommonData &common_data,int tag,bool jacobian = true):
      TetElementForcesAndSourcesCore::UserDataOperator(field_name),
      dAta(data),commonData(common_data),tAg(tag),jAcobian(jacobian) { }

    ublas::vector<adouble> a,dot_W,dp_dt,a_res;
    ublas::matrix<adouble> h,H,invH,F,g;

    vector<double> active;
 
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

	a.resize(3);
	dot_W.resize(3);
	dp_dt.resize(3);
	a_res.resize(3);

	g.resize(3,3);
	h.resize(3,3);
	H.resize(3,3);
	invH.resize(3,3);
	F.resize(3,3);

	dot_W.clear();
	H.clear();
	invH.clear();
	for(int dd = 0;dd<3;dd++) {
	  H(dd,dd) = 1;
	  invH(dd,dd) = 1;
	}

	a_res.resize(3);
	int nb_gauss_pts = row_data.getN().size1();
	commonData.valMass.resize(nb_gauss_pts);
	commonData.jacMassRowPtr.resize(nb_gauss_pts);
	commonData.jacMass.resize(nb_gauss_pts);

	int nb_active_vars = 0;
	for(int gg = 0;gg<nb_gauss_pts;gg++) {

	  if(gg == 0) {

	    trace_on(tAg);
	    
	    for(int nn1 = 0;nn1<3;nn1++) { //0
	      a[nn1] <<= (commonData.dataAtGaussPts["DOT_"+commonData.spatialVelocities][gg])[nn1]; 
	      nb_active_vars++;
	    }
	    for(int nn1 = 0;nn1<3;nn1++) { //3
	      for(int nn2 = 0;nn2<3;nn2++) {
		h(nn1,nn2) <<= (commonData.gradAtGaussPts[commonData.spatialPositions][gg])(nn1,nn2);
		nb_active_vars++;
	      }
	    }
	    if(commonData.dataAtGaussPts["DOT_"+commonData.meshPositions].size()>0) {
	      for(int nn1 = 0;nn1<3;nn1++) { //3+9=12
		for(int nn2 = 0;nn2<3;nn2++) {
		  g(nn1,nn2) <<= (commonData.gradAtGaussPts[commonData.spatialVelocities][gg])(nn1,nn2);
		  nb_active_vars++;
		}
	      }  
	      for(int nn1 = 0;nn1<3;nn1++) { //3+9+9=21
		dot_W(nn1) <<= (commonData.dataAtGaussPts["DOT_"+commonData.meshPositions][gg])[nn1];
		nb_active_vars++;
	      }
	      for(int nn1 = 0;nn1<3;nn1++) { //3+9+9+3=24
		for(int nn2 = 0;nn2<3;nn2++) {
		  H(nn1,nn2) <<= (commonData.gradAtGaussPts[commonData.meshPositions][gg])(nn1,nn2);
		  nb_active_vars++;
		}
	      }
	    }
	    adouble detH = 1;
	    if(commonData.gradAtGaussPts[commonData.meshPositions].size()>0) {
	      ierr = dEterminatnt(H,detH); CHKERRQ(ierr);
	      ierr = iNvert(detH,H,invH); CHKERRQ(ierr);
	    } 

	    double rho0 = dAta.rho0;
	    ublas::vector<double>& a0 = dAta.a0;
	    ierr = calculateFuncUnderIntegral(
	      rho0,a0,a,g,dot_W,H,detH,invH,h,F,dp_dt,a_res); CHKERRQ(ierr);

	    //dependant
	    ublas::vector<double>& res = commonData.valMass[gg];
	    res.resize(3);
	    for(int rr = 0;rr<3;rr++) {
	      a_res[rr] >>= res[rr];
	    }

	    trace_off();

	  }

	  active.resize(nb_active_vars);
	  int aa = 0;
	  for(int nn1 = 0;nn1<3;nn1++) {
	    active[aa++] = (commonData.dataAtGaussPts["DOT_"+commonData.spatialVelocities][gg])[nn1];
	  }
	  for(int nn1 = 0;nn1<3;nn1++) {
	    for(int nn2 = 0;nn2<3;nn2++) {
	      active[aa++] = (commonData.gradAtGaussPts[commonData.spatialPositions][gg])(nn1,nn2);
	    }
	  }
	  if(commonData.dataAtGaussPts["DOT_"+commonData.meshPositions].size()>0) {
	    for(int nn1 = 0;nn1<3;nn1++) {
	      for(int nn2 = 0;nn2<3;nn2++) {
		active[aa++] = (commonData.gradAtGaussPts[commonData.spatialVelocities][gg])(nn1,nn2);
	      }
	    }  
	    for(int nn1 = 0;nn1<3;nn1++) {
	      active[aa++] = (commonData.dataAtGaussPts["DOT_"+commonData.meshPositions][gg])[nn1];
	    }
	    for(int nn1 = 0;nn1<3;nn1++) {
	      for(int nn2 = 0;nn2<3;nn2++) {
		active[aa] = (commonData.gradAtGaussPts[commonData.meshPositions][gg])(nn1,nn2);
	      }
	    }
	  }

	  if(!jAcobian) {
	    ublas::vector<double>& res = commonData.valMass[gg];
	    if(gg>0) {
	      res.resize(3);
	      int r;
	      r = function(tAg,3,nb_active_vars,&active[0],&res[0]);
	      if(r!=3) { // function is locally analytic
		SETERRQ1(PETSC_COMM_SELF,MOFEM_OPERATION_UNSUCCESSFUL,"ADOL-C function evaluation with error r = %d",r);
	      }
	    } 
	    double val = getVolume()*getGaussPts()(3,gg);
	    res *= val;
	  } else {
	    commonData.jacMassRowPtr[gg].resize(3);
	    commonData.jacMass[gg].resize(3,nb_active_vars);
	    for(int nn1 = 0;nn1<3;nn1++) {
	      (commonData.jacMassRowPtr[gg])[nn1] = &(commonData.jacMass[gg](nn1,0));     
	    }
	    int r;
	    r = jacobian(
	      tAg,3,nb_active_vars,
	      &active[0],&(commonData.jacMassRowPtr[gg])[0]);
	    if(r!=3) {
	      SETERRQ(PETSC_COMM_SELF,MOFEM_OPERATION_UNSUCCESSFUL,"ADOL-C function evaluation with error");
	    }
	    double val = getVolume()*getGaussPts()(3,gg);
	    commonData.jacMass[gg] *= val;
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

  struct OpMassRhs: public TetElementForcesAndSourcesCore::UserDataOperator,CommonFunctions {

    BlockData &dAta;
    CommonData &commonData;

    OpMassRhs(const string field_name,BlockData &data,CommonData &common_data):
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
	  ublas::vector<double>& res = commonData.valMass[gg];
	  //cerr << res << endl;
	  for(int dd = 0;dd<nb_dofs/3;dd++) {
	    for(int rr = 0;rr<3;rr++) {
	      nf[3*dd+rr] += row_data.getN()(gg,dd)*res[rr];
	    }
	  }
	}

	if((unsigned int)nb_dofs > 3*row_data.getN().size2()) {
	  SETERRQ(PETSC_COMM_SELF,1,"data inconsistency");
	}
	ierr = VecSetValues(getFEMethod()->ts_F,nb_dofs,
	  &row_data.getIndices()[0],&nf[0],ADD_VALUES); CHKERRQ(ierr);

      } catch (const std::exception& ex) {
	ostringstream ss;
	ss << "throw in method: " << ex.what() << endl;
	SETERRQ(PETSC_COMM_SELF,1,ss.str().c_str());
      }

      PetscFunctionReturn(0);
    }

  };

  struct OpMassLhs_dM_dv: public TetElementForcesAndSourcesCore::UserDataOperator,CommonFunctions {

    BlockData &dAta;
    CommonData &commonData;
    int tAg;

    OpMassLhs_dM_dv(const string vel_field,const string field_name,BlockData &data,CommonData &common_data):
      TetElementForcesAndSourcesCore::UserDataOperator(vel_field,field_name),
      dAta(data),commonData(common_data) { symm = false;  }

    ublas::matrix<double> k,jac;

    virtual PetscErrorCode getJac(DataForcesAndSurcesCore::EntData &col_data,int gg) {
      PetscFunctionBegin;
      int nb_col = col_data.getIndices().size();
      jac.clear();
      //cerr << commonData.jacMass[gg] << endl;
      //cerr << jac << endl;
      ublas::vector<double> N = col_data.getN(gg,nb_col/3);
      for(int dd = 0;dd<nb_col/3;dd++) {
	for(int nn = 0;nn<3;nn++) {
	  jac(0,3*dd+nn) = commonData.jacMass[gg](0,nn)*N(dd)*getFEMethod()->ts_a; 
	  jac(1,3*dd+nn) = commonData.jacMass[gg](1,nn)*N(dd)*getFEMethod()->ts_a; 
	  jac(2,3*dd+nn) = commonData.jacMass[gg](2,nn)*N(dd)*getFEMethod()->ts_a; 
	}
      }
      if(commonData.dataAtGaussPts["DOT_"+commonData.meshPositions].size()>0) {
	ublas::matrix<double> diffN = col_data.getDiffN(gg,nb_col/3);
	for(int dd = 0;dd<nb_col/3;dd++) {
	  //h00 //h01 //h02
	  jac(0,3*dd+0) += commonData.jacMass[gg](0,3+9+3*0+0)*diffN(dd,0);
	  jac(0,3*dd+0) += commonData.jacMass[gg](0,3+9+3*0+1)*diffN(dd,1);
	  jac(0,3*dd+0) += commonData.jacMass[gg](0,3+9+3*0+2)*diffN(dd,2);
	  //h10 //h11 //h12
	  jac(1,3*dd+1) += commonData.jacMass[gg](1,3+9+3*1+0)*diffN(dd,0);
	  jac(1,3*dd+1) += commonData.jacMass[gg](1,3+9+3*1+1)*diffN(dd,1);
	  jac(1,3*dd+1) += commonData.jacMass[gg](1,3+9+3*1+2)*diffN(dd,2);
	  //h20 //h21 //h22
	  jac(2,3*dd+2) += commonData.jacMass[gg](2,3+9+3*2+0)*diffN(dd,0);
	  jac(2,3*dd+2) += commonData.jacMass[gg](2,3+9+3*2+1)*diffN(dd,1);
	  jac(2,3*dd+2) += commonData.jacMass[gg](2,3+9+3*2+2)*diffN(dd,2);
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
      if(nb_row==0) PetscFunctionReturn(0);
      if(nb_col==0) PetscFunctionReturn(0);

      try {

	k.resize(nb_row,nb_col);
	k.clear();
	jac.resize(3,nb_col);

	for(unsigned int gg = 0;gg<row_data.getN().size1();gg++) {

	  ierr = getJac(col_data,gg); CHKERRQ(ierr);

	  { //integrate element stiffnes matrix
	    for(int dd1 = 0;dd1<nb_row/3;dd1++) {
	      for(int rr1 = 0;rr1<3;rr1++) {
		for(int dd2 = 0;dd2<nb_col/3;dd2++) {
		  for(int rr2 = 0;rr2<3;rr2++) {
		    k(3*dd1+rr1,3*dd2+rr2) += row_data.getN()(gg,dd1)*jac(rr1,3*dd2+rr2);
		  }
		}
	      }
	    }
	  }

	}

	ierr = MatSetValues(getFEMethod()->ts_B,
	  nb_row,&row_data.getIndices()[0],
	  nb_col,&col_data.getIndices()[0],
	  &k(0,0),ADD_VALUES); CHKERRQ(ierr);

      } catch (const std::exception& ex) {
	ostringstream ss;
	ss << "throw in method: " << ex.what() << endl;
	SETERRQ(PETSC_COMM_SELF,1,ss.str().c_str());
      }

      PetscFunctionReturn(0);
    }

  };

  struct OpMassLhs_dM_dx: public OpMassLhs_dM_dv  {

    OpMassLhs_dM_dx(const string field_name,const string col_field,BlockData &data,CommonData &common_data):
      OpMassLhs_dM_dv(field_name,col_field,data,common_data) {}

    PetscErrorCode getJac(DataForcesAndSurcesCore::EntData &col_data,int gg) {
      PetscFunctionBegin;
      int nb_col = col_data.getIndices().size();
      jac.clear();
      ublas::matrix<double> diffN = col_data.getDiffN(gg,nb_col/3);
      for(int dd = 0;dd<nb_col/3;dd++) {
	//h00 //h01 //h02
	jac(0,3*dd+0) += commonData.jacMass[gg](0,3+3*0+0)*diffN(dd,0);
	jac(0,3*dd+0) += commonData.jacMass[gg](0,3+3*0+1)*diffN(dd,1);
	jac(0,3*dd+0) += commonData.jacMass[gg](0,3+3*0+2)*diffN(dd,2);
	//h10 //h11 //h12
	jac(1,3*dd+1) += commonData.jacMass[gg](1,3+3*1+0)*diffN(dd,0);
	jac(1,3*dd+1) += commonData.jacMass[gg](1,3+3*1+1)*diffN(dd,1);
	jac(1,3*dd+1) += commonData.jacMass[gg](1,3+3*1+2)*diffN(dd,2);
	//h20 //h21 //h22
	jac(2,3*dd+2) += commonData.jacMass[gg](2,3+3*2+0)*diffN(dd,0);
	jac(2,3*dd+2) += commonData.jacMass[gg](2,3+3*2+1)*diffN(dd,1);
	jac(2,3*dd+2) += commonData.jacMass[gg](2,3+3*2+2)*diffN(dd,2);
      }
      PetscFunctionReturn(0);
    }

  };

  struct OpMassLhs_dM_dX: public OpMassLhs_dM_dv  {

    OpMassLhs_dM_dX(const string field_name,const string col_field,BlockData &data,CommonData &common_data):
      OpMassLhs_dM_dv(field_name,col_field,data,common_data) {}

    PetscErrorCode getJac(DataForcesAndSurcesCore::EntData &col_data,int gg) {
      PetscFunctionBegin;
      int nb_col = col_data.getIndices().size();
      jac.clear();
      //cerr << commonData.jacVel[gg] << endl;
      //cerr << jac << endl;
      ublas::vector<double> N = col_data.getN(gg,nb_col/3);
      for(int dd = 0;dd<nb_col/3;dd++) {
	for(int nn = 0;nn<3;nn++) {
	  jac(0,3*dd+nn) = commonData.jacVel[gg](0,3+9+9+nn)*N(dd)*getFEMethod()->ts_a; 
	  jac(1,3*dd+nn) = commonData.jacVel[gg](1,3+9+9+nn)*N(dd)*getFEMethod()->ts_a; 
	  jac(2,3*dd+nn) = commonData.jacVel[gg](2,3+9+9+nn)*N(dd)*getFEMethod()->ts_a; 
	}
      }
      ublas::matrix<double> diffN = col_data.getDiffN(gg,nb_col/3);
      for(int dd = 0;dd<nb_col/3;dd++) {
	//h00 //h01 //h02
	jac(0,3*dd+0) += commonData.jacVel[gg](0,3+9+9+3+3*0+0)*diffN(dd,0);
	jac(0,3*dd+0) += commonData.jacVel[gg](0,3+9+9+3+3*0+1)*diffN(dd,1);
	jac(0,3*dd+0) += commonData.jacVel[gg](0,3+9+9+3+3*0+2)*diffN(dd,2);
	//h10 //h11 //h12
	jac(1,3*dd+1) += commonData.jacVel[gg](1,3+9+9+3+3*1+0)*diffN(dd,0);
	jac(1,3*dd+1) += commonData.jacVel[gg](1,3+9+9+3+3*1+1)*diffN(dd,1);
	jac(1,3*dd+1) += commonData.jacVel[gg](1,3+9+9+3+3*1+2)*diffN(dd,2);
	//h20 //h21 //h22
	jac(2,3*dd+2) += commonData.jacVel[gg](2,3+9+9+3+3*2+0)*diffN(dd,0);
	jac(2,3*dd+2) += commonData.jacVel[gg](2,3+9+9+3+3*2+1)*diffN(dd,1);
	jac(2,3*dd+2) += commonData.jacVel[gg](2,3+9+9+3+3*2+2)*diffN(dd,2);
      }
      PetscFunctionReturn(0);
    }

  };

  struct OpVelocityJacobian: public TetElementForcesAndSourcesCore::UserDataOperator,CommonFunctions {

    BlockData &dAta;
    CommonData &commonData;
    int tAg;
    bool jAcobian;

    OpVelocityJacobian(const string field_name,BlockData &data,CommonData &common_data,int tag,bool jacobian = true):
      TetElementForcesAndSourcesCore::UserDataOperator(field_name),
      dAta(data),commonData(common_data),tAg(tag),jAcobian(jacobian) { }

    ublas::vector<adouble> a_res;
    ublas::vector<adouble> v,dot_w,dot_W;
    ublas::matrix<adouble> h,H,invH,F;
    ublas::vector<adouble> dot_u;
    adouble detH;

    vector<double> active;
 
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

	v.resize(3);
	dot_w.resize(3);
	h.resize(3,3);
	h.clear();
	F.resize(3,3);
	dot_W.resize(3);
	dot_W.clear();
	H.resize(3,3);
	H.clear();
	invH.resize(3,3);
	invH.clear();
	dot_u.resize(3);
	for(int dd = 0;dd<3;dd++) {
	  H(dd,dd) = 1;
	  invH(dd,dd) = 1;     
	}

	a_res.resize(3);
	int nb_gauss_pts = row_data.getN().size1();
	commonData.valVel.resize(nb_gauss_pts);
	commonData.jacVelRowPtr.resize(nb_gauss_pts);
	commonData.jacVel.resize(nb_gauss_pts);

	int nb_active_vars = 0;
	for(int gg = 0;gg<nb_gauss_pts;gg++) {

	  if(gg == 0) {

	    trace_on(tAg);

	    for(int nn1 = 0;nn1<3;nn1++) { //0
	      v[nn1] <<= commonData.dataAtGaussPts[commonData.spatialVelocities][gg][nn1]; nb_active_vars++;
	    }
	    for(int nn1 = 0;nn1<3;nn1++) { //3
	      dot_w[nn1] <<= commonData.dataAtGaussPts["DOT_"+commonData.spatialPositions][gg][nn1]; nb_active_vars++;
	    }
	    if(commonData.dataAtGaussPts["DOT_"+commonData.meshPositions].size()>0) {
	      for(int nn1 = 0;nn1<3;nn1++) { //3+3 = 6
		for(int nn2 = 0;nn2<3;nn2++) {
		  h(nn1,nn2) <<= commonData.gradAtGaussPts[commonData.spatialPositions][gg](nn1,nn2); nb_active_vars++;
		}
	      }
	      for(int nn1 = 0;nn1<3;nn1++) { //3+3+9
		dot_W[nn1] <<= commonData.dataAtGaussPts["DOT_"+commonData.meshPositions][gg][nn1]; nb_active_vars++;
	      }
	    }
	    if(commonData.gradAtGaussPts[commonData.meshPositions].size()>0) {
	      for(int nn1 = 0;nn1<3;nn1++) { //3+3+9+3
		for(int nn2 = 0;nn2<3;nn2++) {
		  H(nn1,nn2) <<= commonData.gradAtGaussPts[commonData.meshPositions][gg](nn1,nn2); nb_active_vars++;
		}
	      }
	    }
	    detH = 1;
	    if(commonData.gradAtGaussPts[commonData.meshPositions].size()>0) {
	      ierr = dEterminatnt(H,detH); CHKERRQ(ierr);
	    }
	    ierr = iNvert(detH,H,invH); CHKERRQ(ierr);
	    noalias(F) = prod(h,invH);
	    noalias(dot_u) = dot_w - prod(F,dot_W);
	    noalias(a_res) = (v - dot_u)*detH;
	    //dependant
	    ublas::vector<double>& res = commonData.valVel[gg];
	    res.resize(3);
	    for(int rr = 0;rr<3;rr++) {
	      a_res[rr] >>= res[rr];
	    }
	    trace_off();

	  }

	  active.resize(nb_active_vars);
	  int aa = 0;
	  for(int nn1 = 0;nn1<3;nn1++) {
	    active[aa++] = commonData.dataAtGaussPts[commonData.spatialVelocities][gg][nn1]; 
	  }
	  for(int nn1 = 0;nn1<3;nn1++) {
	    active[aa++] = commonData.dataAtGaussPts["DOT_"+commonData.spatialPositions][gg][nn1]; 
	  }
	  if(commonData.dataAtGaussPts["DOT_"+commonData.meshPositions].size()>0) {
	    for(int nn1 = 0;nn1<3;nn1++) {
	      for(int nn2 = 0;nn2<3;nn2++) {
		active[aa++] = commonData.gradAtGaussPts[commonData.spatialPositions][gg](nn1,nn2);
	      }	
	    }
	    for(int nn1 = 0;nn1<3;nn1++) {
	      active[aa++] = commonData.dataAtGaussPts["DOT_"+commonData.meshPositions][gg][nn1];
	    }
	  }
	  if(commonData.gradAtGaussPts[commonData.meshPositions].size()>0) {
	    for(int nn1 = 0;nn1<3;nn1++) {
	      for(int nn2 = 0;nn2<3;nn2++) {
		active[aa++] = commonData.gradAtGaussPts[commonData.meshPositions][gg](nn1,nn2);
	      }
	    }
	  }

	  if(!jAcobian) {
	    ublas::vector<double>& res = commonData.valVel[gg];
	    if(gg>0) {
	      res.resize(3);
	      int r;
	      r = function(tAg,3,nb_active_vars,&active[0],&res[0]);
	      if(r!=3) {
		SETERRQ(PETSC_COMM_SELF,MOFEM_OPERATION_UNSUCCESSFUL,"ADOL-C function evaluation with error");
	      }
	    } 
	    double val = getVolume()*getGaussPts()(3,gg);
	    res *= val;
	  } else {
	    commonData.jacVelRowPtr[gg].resize(3);
	    commonData.jacVel[gg].resize(3,nb_active_vars);
	    for(int nn1 = 0;nn1<3;nn1++) {
	      (commonData.jacVelRowPtr[gg])[nn1] = &(commonData.jacVel[gg](nn1,0));     
	    }
	    int r;
	    r = jacobian(
	      tAg,3,nb_active_vars,
	      &active[0],&(commonData.jacVelRowPtr[gg])[0]);
	    if(r!=3) {
	      SETERRQ(PETSC_COMM_SELF,MOFEM_OPERATION_UNSUCCESSFUL,"ADOL-C function evaluation with error");
	    }
	    double val = getVolume()*getGaussPts()(3,gg);
	    commonData.jacVel[gg] *= val;
	    //cerr << gg << " : " << commonData.jacVel[gg] << endl;
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

  struct OpVelocityRhs: public TetElementForcesAndSourcesCore::UserDataOperator,CommonFunctions {

    BlockData &dAta;
    CommonData &commonData;

    OpVelocityRhs(const string field_name,BlockData &data,CommonData &common_data):
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
      int nb_dofs = row_data.getIndices().size();
      if(nb_dofs==0) PetscFunctionReturn(0);

      try {

	nf.resize(nb_dofs);
	nf.clear();

	for(unsigned int gg = 0;gg<row_data.getN().size1();gg++) {
	  ublas::vector<double>& res = commonData.valVel[gg];
	  for(int dd = 0;dd<nb_dofs/3;dd++) {
	    for(int rr = 0;rr<3;rr++) {
	      nf[3*dd+rr] += row_data.getN()(gg,dd)*res[rr];
	    }
	  }
	}

	if(row_data.getIndices().size() > 3*row_data.getN().size2()) {
	  SETERRQ(PETSC_COMM_SELF,1,"data inconsistency");
	}
	ierr = VecSetValues(getFEMethod()->ts_F,row_data.getIndices().size(),
	  &row_data.getIndices()[0],&nf[0],ADD_VALUES); CHKERRQ(ierr);

      } catch (const std::exception& ex) {
	ostringstream ss;
	ss << "throw in method: " << ex.what() << endl;
	SETERRQ(PETSC_COMM_SELF,1,ss.str().c_str());
      }

      PetscFunctionReturn(0);
    }

  };

  struct OpVelocityLhs_dV_dv: public OpMassLhs_dM_dv {

    OpVelocityLhs_dV_dv(const string vel_field,const string field_name,BlockData &data,CommonData &common_data):
      OpMassLhs_dM_dv(vel_field,field_name,data,common_data) {};

    virtual PetscErrorCode getJac(DataForcesAndSurcesCore::EntData &col_data,int gg) {
      PetscFunctionBegin;
      int nb_col = col_data.getIndices().size();
      jac.clear();
      ublas::vector<double> N = col_data.getN(gg,nb_col/3);
      //cerr << commonData.jacVel[gg] << endl;
      for(int dd = 0;dd<nb_col/3;dd++) {
	for(int nn = 0;nn<3;nn++) {
	  jac(0,3*dd+nn) = commonData.jacVel[gg](0,nn)*N(dd); 
	  jac(1,3*dd+nn) = commonData.jacVel[gg](1,nn)*N(dd); 
	  jac(2,3*dd+nn) = commonData.jacVel[gg](2,nn)*N(dd); 
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
      if(nb_row==0) PetscFunctionReturn(0);
      if(nb_col==0) PetscFunctionReturn(0);

      try {

	k.resize(nb_row,nb_col);
	k.clear();
	jac.resize(3,nb_col);

	for(unsigned int gg = 0;gg<row_data.getN().size1();gg++) {

	  ierr = getJac(col_data,gg); CHKERRQ(ierr);

	  { //integrate element stiffnes matrix
	    for(int dd1 = 0;dd1<nb_row/3;dd1++) {
	      for(int rr1 = 0;rr1<3;rr1++) {
		for(int dd2 = 0;dd2<nb_col/3;dd2++) {
		  for(int rr2 = 0;rr2<3;rr2++) {
		    k(3*dd1+rr1,3*dd2+rr2) += row_data.getN()(gg,dd1)*jac(rr1,3*dd2+rr2);
		  }
		}
	      }
	    }
	  }

	}

	ierr = MatSetValues(getFEMethod()->ts_B,
	  nb_row,&row_data.getIndices()[0],
	  nb_col,&col_data.getIndices()[0],
	  &k(0,0),ADD_VALUES); CHKERRQ(ierr);


      } catch (const std::exception& ex) {
	ostringstream ss;
	ss << "throw in method: " << ex.what() << endl;
	SETERRQ(PETSC_COMM_SELF,1,ss.str().c_str());
      }

      PetscFunctionReturn(0);
    }

  };

  struct OpVelocityLhs_dV_dx: public OpVelocityLhs_dV_dv {

    OpVelocityLhs_dV_dx(const string vel_field,const string field_name,BlockData &data,CommonData &common_data):
      OpVelocityLhs_dV_dv(vel_field,field_name,data,common_data) {}

    virtual PetscErrorCode getJac(DataForcesAndSurcesCore::EntData &col_data,int gg) {
      PetscFunctionBegin;
      int nb_col = col_data.getIndices().size();
      jac.clear();
      //cerr << commonData.jacVel[gg] << endl;
      //cerr << jac << endl;
      ublas::vector<double> N = col_data.getN(gg,nb_col/3);
      for(int dd = 0;dd<nb_col/3;dd++) {
	for(int nn = 0;nn<3;nn++) {
	  jac(0,3*dd+nn) = commonData.jacVel[gg](0,3+nn)*N(dd)*getFEMethod()->ts_a; 
	  jac(1,3*dd+nn) = commonData.jacVel[gg](1,3+nn)*N(dd)*getFEMethod()->ts_a; 
	  jac(2,3*dd+nn) = commonData.jacVel[gg](2,3+nn)*N(dd)*getFEMethod()->ts_a; 
	}
      }
      if(commonData.dataAtGaussPts["DOT_"+commonData.meshPositions].size()>0) {
        ublas::matrix<double> diffN = col_data.getDiffN(gg,nb_col/3);
        for(int dd = 0;dd<nb_col/3;dd++) {
	  //h00 //h01 //h02
	  jac(0,3*dd+0) += commonData.jacVel[gg](0,3+3+3*0+0)*diffN(dd,0);
	  jac(0,3*dd+0) += commonData.jacVel[gg](0,3+3+3*0+1)*diffN(dd,1);
	  jac(0,3*dd+0) += commonData.jacVel[gg](0,3+3+3*0+2)*diffN(dd,2);
	  //h10 //h11 //h12
	  jac(1,3*dd+1) += commonData.jacVel[gg](1,3+3+3*1+0)*diffN(dd,0);
	  jac(1,3*dd+1) += commonData.jacVel[gg](1,3+3+3*1+1)*diffN(dd,1);
	  jac(1,3*dd+1) += commonData.jacVel[gg](1,3+3+3*1+2)*diffN(dd,2);
	  //h20 //h21 //h22
	  jac(2,3*dd+2) += commonData.jacVel[gg](2,3+3+3*2+0)*diffN(dd,0);
	  jac(2,3*dd+2) += commonData.jacVel[gg](2,3+3+3*2+1)*diffN(dd,1);
	  jac(2,3*dd+2) += commonData.jacVel[gg](2,3+3+3*2+2)*diffN(dd,2);
        }
      }
      //cerr << row_field_name << " " << col_field_name << endl;
      PetscFunctionReturn(0);
    }

  };


  struct OpVelocityLhs_dV_dX: public OpVelocityLhs_dV_dv {

    OpVelocityLhs_dV_dX(const string vel_field,const string field_name,BlockData &data,CommonData &common_data):
      OpVelocityLhs_dV_dv(vel_field,field_name,data,common_data) {}

    virtual PetscErrorCode getJac(DataForcesAndSurcesCore::EntData &col_data,int gg) {
      PetscFunctionBegin;
      int nb_col = col_data.getIndices().size();
      jac.clear();
      //cerr << commonData.jacVel[gg] << endl;
      //cerr << jac << endl;
      ublas::vector<double> N = col_data.getN(gg,nb_col/3);
      for(int dd = 0;dd<nb_col/3;dd++) {
	for(int nn = 0;nn<3;nn++) {
	  jac(0,3*dd+nn) = commonData.jacVel[gg](0,3+3+9+nn)*N(dd)*getFEMethod()->ts_a; 
	  jac(1,3*dd+nn) = commonData.jacVel[gg](1,3+3+9+nn)*N(dd)*getFEMethod()->ts_a; 
	  jac(2,3*dd+nn) = commonData.jacVel[gg](2,3+3+9+nn)*N(dd)*getFEMethod()->ts_a; 
	}
      }
      ublas::matrix<double> diffN = col_data.getDiffN(gg,nb_col/3);
      for(int dd = 0;dd<nb_col/3;dd++) {
	//h00 //h01 //h02
	jac(0,3*dd+0) += commonData.jacVel[gg](0,3+3+9+3+3*0+0)*diffN(dd,0);
	jac(0,3*dd+0) += commonData.jacVel[gg](0,3+3+9+3+3*0+1)*diffN(dd,1);
	jac(0,3*dd+0) += commonData.jacVel[gg](0,3+3+9+3+3*0+2)*diffN(dd,2);
	//h10 //h11 //h12
	jac(1,3*dd+1) += commonData.jacVel[gg](1,3+3+9+3+3*1+0)*diffN(dd,0);
	jac(1,3*dd+1) += commonData.jacVel[gg](1,3+3+9+3+3*1+1)*diffN(dd,1);
	jac(1,3*dd+1) += commonData.jacVel[gg](1,3+3+9+3+3*1+2)*diffN(dd,2);
	//h20 //h21 //h22
	jac(2,3*dd+2) += commonData.jacVel[gg](2,3+3+9+3+3*2+0)*diffN(dd,0);
	jac(2,3*dd+2) += commonData.jacVel[gg](2,3+3+9+3+3*2+1)*diffN(dd,1);
	jac(2,3*dd+2) += commonData.jacVel[gg](2,3+3+9+3+3*2+2)*diffN(dd,2);
      }

      //cerr << row_field_name << " " << col_field_name << endl;

      PetscFunctionReturn(0);
    }

  };

  struct OpEshelbyDynamicMaterialMomentumJacobian: public TetElementForcesAndSourcesCore::UserDataOperator,CommonFunctions {

    BlockData &dAta;
    CommonData &commonData;
    int tAg;
    bool jAcobian;

    OpEshelbyDynamicMaterialMomentumJacobian(const string field_name,BlockData &data,CommonData &common_data,int tag,bool jacobian = true):
      TetElementForcesAndSourcesCore::UserDataOperator(field_name),
      dAta(data),commonData(common_data),tAg(tag),jAcobian(jacobian) {}

    ublas::vector<adouble> v;
    ublas::matrix<adouble> H,invH,h,F;
    ublas::vector<double> active;

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

	v.resize(3);
	h.resize(3,3);	
	F.resize(3,3);
	H.resize(3,3);
	H.clear();
	invH.resize(3,3);
	invH.clear();
	for(int dd = 0;dd<3;dd++) {
	  H(dd,dd) = 1;
	  invH(dd,dd) = 1;     
	}

	int nb_gauss_pts = row_data.getN().size1();
	commonData.valT.resize(nb_gauss_pts);
	commonData.jacEshelby.resize(nb_gauss_pts);

	int nb_active_vars = 0;
	for(int gg = 0;gg<nb_gauss_pts;gg++) {

	  if(gg == 0) {

	    trace_on(tAg);

	    for(int nn1 = 0;nn1<3;nn1++) { //0
	      v[nn1] <<= commonData.dataAtGaussPts[commonData.spatialVelocities][gg][nn1]; nb_active_vars++;
	    }
	    for(int nn1 = 0;nn1<3;nn1++) { //3
	      for(int nn2 = 0;nn2<3;nn2++) {
		h(nn1,nn2) <<= commonData.gradAtGaussPts[commonData.spatialPositions][gg](nn1,nn2); nb_active_vars++;
	      }
	    }
	    if(commonData.gradAtGaussPts[commonData.meshPositions].size()>0) {
	      for(int nn1 = 0;nn1<3;nn1++) { //3+9
		for(int nn2 = 0;nn2<3;nn2++) {
		  H(nn1,nn2) <<= commonData.gradAtGaussPts[commonData.meshPositions][gg](nn1,nn2); nb_active_vars++;
		}
	      }
	    }
	    adouble detH;
	    detH = 1;
	    if(commonData.gradAtGaussPts[commonData.meshPositions].size()>0) {
	      ierr = dEterminatnt(H,detH); CHKERRQ(ierr);
	    }
	    ierr = iNvert(detH,H,invH); CHKERRQ(ierr);
	    noalias(F) = prod(h,invH);
	    adouble detF;
	    ierr = dEterminatnt(F,detF); CHKERRQ(ierr);
	    double rho0 = dAta.rho0;
	    adouble rho = rho0*detF;
	    adouble a_T = 0.5*rho*inner_prod(v,v)*detH;
	    //dependant
	    a_T >>= commonData.valT[gg];
	    trace_off();

	  }

	  active.resize(nb_active_vars);
	  int aa = 0;
	  for(int nn1 = 0;nn1<3;nn1++) {
	    active[aa++] = commonData.dataAtGaussPts[commonData.spatialVelocities][gg][nn1]; 
	  }
	  for(int nn1 = 0;nn1<3;nn1++) { 
	    for(int nn2 = 0;nn2<3;nn2++) {
	      active[aa++] = commonData.gradAtGaussPts[commonData.spatialPositions][gg](nn1,nn2); 
	    }
	  }
	  if(commonData.gradAtGaussPts[commonData.meshPositions].size()>0) {
	    for(int nn1 = 0;nn1<3;nn1++) {
	      for(int nn2 = 0;nn2<3;nn2++) {
		active[aa++] = commonData.gradAtGaussPts[commonData.meshPositions][gg](nn1,nn2);
	      }
	    }
	  }

	  if(!jAcobian) {
	    if(gg>0) {
	      int r;
	      r = function(tAg,1,nb_active_vars,&active[0],&commonData.valT[gg]);
	      if(r!=3) {
		SETERRQ(PETSC_COMM_SELF,MOFEM_OPERATION_UNSUCCESSFUL,"ADOL-C function evaluation with error");
	      }
	    } 
	    double val = getVolume()*getGaussPts()(3,gg);
	    commonData.valT[gg] *= val;
	  } else {
	    commonData.jacEshelby[gg].resize(nb_active_vars);
	    commonData.jacEshelby[gg].clear();
	    int r;
	    r = gradient(
	      tAg,nb_active_vars,
	      &active[0],&(commonData.jacEshelby[gg][0]));
	    if(r!=3) {
	      SETERRQ(PETSC_COMM_SELF,MOFEM_OPERATION_UNSUCCESSFUL,"ADOL-C function evaluation with error");
	    }
	    double val = getVolume()*getGaussPts()(3,gg);
	    commonData.jacEshelby[gg] *= val;
	    //cerr << gg << " : " << commonData.jacVel[gg] << endl;
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

  struct OpEshelbyDynamicMaterialMomentumRhs: public TetElementForcesAndSourcesCore::UserDataOperator,CommonFunctions {

    BlockData &dAta;
    CommonData &commonData;

    OpEshelbyDynamicMaterialMomentumRhs(const string field_name,BlockData &data,CommonData &common_data):
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
      int nb_dofs = row_data.getIndices().size();
      if(nb_dofs==0) PetscFunctionReturn(0);

      try {

	nf.resize(nb_dofs);
	nf.clear();

	for(unsigned int gg = 0;gg<row_data.getN().size1();gg++) {
	  const DataForcesAndSurcesCore::MatrixAdaptor &diffN = row_data.getDiffN(gg,nb_dofs/3);
	  double mT = -commonData.valT[gg];
	  for(int dd = 0;dd<nb_dofs/3;dd++) {
	    for(int rr = 0;rr<3;rr++) {
	      nf[3*dd+rr] += diffN(dd,rr)*mT;
	    }
	  }
	}

	if(row_data.getIndices().size() > 3*row_data.getN().size2()) {
	  SETERRQ(PETSC_COMM_SELF,1,"data inconsistency");
	}
	ierr = VecSetValues(getFEMethod()->ts_F,row_data.getIndices().size(),
	  &row_data.getIndices()[0],&nf[0],ADD_VALUES); CHKERRQ(ierr);

      } catch (const std::exception& ex) {
	ostringstream ss;
	ss << "throw in method: " << ex.what() << endl;
	SETERRQ(PETSC_COMM_SELF,1,ss.str().c_str());
      }

      PetscFunctionReturn(0);
    }

  };

  struct OpEshelbyDynamicMaterialMomentumRhs_dv: public OpMassLhs_dM_dv {

    OpEshelbyDynamicMaterialMomentumRhs_dv(const string vel_field,const string field_name,BlockData &data,CommonData &common_data):
      OpMassLhs_dM_dv(vel_field,field_name,data,common_data) {};

    ublas::vector<double> jac;

    virtual PetscErrorCode getJac(DataForcesAndSurcesCore::EntData &col_data,int gg) {
      PetscFunctionBegin;
      int nb_col = col_data.getIndices().size();
      jac.clear();
      ublas::vector<double> N = col_data.getN(gg,nb_col/3);
      //cerr << commonData.jacVel[gg] << endl;
      for(int dd = 0;dd<nb_col/3;dd++) {
	for(int nn = 0;nn<3;nn++) {
	  jac(3*dd+nn) = commonData.jacEshelby[gg](nn)*N(dd); 
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
      if(nb_row==0) PetscFunctionReturn(0);
      if(nb_col==0) PetscFunctionReturn(0);

      try {

	k.resize(nb_row,nb_col);
	k.clear();
	jac.resize(nb_col);

	for(unsigned int gg = 0;gg<row_data.getN().size1();gg++) {

	  ierr = getJac(col_data,gg); CHKERRQ(ierr);

	  const DataForcesAndSurcesCore::MatrixAdaptor &diffN = row_data.getDiffN(gg,nb_row/3);

	  { //integrate element stiffnes matrix
	    for(int dd1 = 0;dd1<nb_row/3;dd1++) {
	      for(int rr1 = 0;rr1<3;rr1++) {
		for(int dd2 = 0;dd2<nb_col/3;dd2++) {
		  for(int rr2 = 0;rr2<3;rr2++) {
		    k(3*dd1+rr1,3*dd2+rr2) += diffN(dd1,rr2)*jac(3*dd2+rr2);
		  }
		}
	      }
	    }
	  }

	}

	ierr = MatSetValues(getFEMethod()->ts_B,
	  nb_row,&row_data.getIndices()[0],
	  nb_col,&col_data.getIndices()[0],
	  &k(0,0),ADD_VALUES); CHKERRQ(ierr);


      } catch (const std::exception& ex) {
	ostringstream ss;
	ss << "throw in method: " << ex.what() << endl;
	SETERRQ(PETSC_COMM_SELF,1,ss.str().c_str());
      }

      PetscFunctionReturn(0);
    }

  };

  struct OpEshelbyDynamicMaterialMomentumRhs_dx: public OpEshelbyDynamicMaterialMomentumRhs_dv {

    OpEshelbyDynamicMaterialMomentumRhs_dx(const string vel_field,const string field_name,BlockData &data,CommonData &common_data):
      OpEshelbyDynamicMaterialMomentumRhs_dv(vel_field,field_name,data,common_data) {};

    ublas::vector<double> jac;

    virtual PetscErrorCode getJac(DataForcesAndSurcesCore::EntData &col_data,int gg) {
      PetscFunctionBegin;
      jac.clear();
      int nb_col = col_data.getFieldData().size();
      const DataForcesAndSurcesCore::MatrixAdaptor diffN = col_data.getDiffN(gg,nb_col/3);
      for(int dd = 0;dd<nb_col/3;dd++) {
	for(int rr = 0;rr<3;rr++) {
	  for(int jj = 0;jj<3;jj++) {
	    jac(3*dd+rr) += commonData.jacEshelby[gg](3+3*rr+jj)*diffN(dd,jj);
	  }
	}	
      }
      PetscFunctionReturn(0);
    }

  };

  struct OpEshelbyDynamicMaterialMomentumRhs_dX: public OpEshelbyDynamicMaterialMomentumRhs_dv {

    OpEshelbyDynamicMaterialMomentumRhs_dX(const string vel_field,const string field_name,BlockData &data,CommonData &common_data):
      OpEshelbyDynamicMaterialMomentumRhs_dv(vel_field,field_name,data,common_data) {};

    ublas::vector<double> jac;

    virtual PetscErrorCode getJac(DataForcesAndSurcesCore::EntData &col_data,int gg) {
      PetscFunctionBegin;
      jac.clear();
      int nb_col = col_data.getFieldData().size();
      const DataForcesAndSurcesCore::MatrixAdaptor diffN = col_data.getDiffN(gg,nb_col/3);
      for(int dd = 0;dd<nb_col/3;dd++) {
	for(int rr = 0;rr<3;rr++) {
	  for(int jj = 0;jj<3;jj++) {
	    jac(3*dd+rr) += commonData.jacEshelby[gg](3+9+3*rr+jj)*diffN(dd,jj);
	  }
	}	
      }
      PetscFunctionReturn(0);
    }

  };

  struct UpdateAndControl: public FEMethod {

    FieldInterface& mField;
    TS tS;
    const string velocityField;
    const string spatialPositionField;

    int jacobianLag;
    UpdateAndControl(FieldInterface& _mField,TS _ts,
      const string velocity_field,
      const string spatial_position_field):
      mField(_mField),tS(_ts),
      velocityField(velocity_field),
      spatialPositionField(spatial_position_field),
      jacobianLag(-1) {}

    PetscErrorCode preProcess() {
      PetscFunctionBegin;
      PetscErrorCode ierr;

      switch (ts_ctx) {
	case CTX_TSSETIFUNCTION: {
	  snes_ctx = CTX_SNESSETFUNCTION;
	  snes_f = ts_F;
	  break;
	}
	case CTX_TSSETIJACOBIAN: {
	  snes_ctx = CTX_SNESSETJACOBIAN;
	  snes_B = ts_B;
	  break;
	}
	default:
	break;
      }

      ierr = mField.set_other_local_VecCreateGhost(problemPtr,velocityField,"DOT_"+velocityField,COL,ts_u_t,INSERT_VALUES,SCATTER_REVERSE); CHKERRQ(ierr);
      ierr = mField.set_other_local_VecCreateGhost(problemPtr,spatialPositionField,"DOT_"+spatialPositionField,COL,ts_u_t,INSERT_VALUES,SCATTER_REVERSE); CHKERRQ(ierr);
      PetscFunctionReturn(0);
    }

    PetscErrorCode postProcess() {
      PetscFunctionBegin;
      //PetscErrorCode ierr;
      //SNES snes;
      //ierr = TSGetSNES(tS,&snes); CHKERRQ(ierr);
      //ierr = SNESSetLagJacobian(snes,jacobianLag); CHKERRQ(ierr);
      PetscFunctionReturn(0);
    }

  };


  PetscErrorCode setBlocks(bool get_density_form_elastic_block_set = true) {
    PetscFunctionBegin;
    ErrorCode rval;
    PetscErrorCode ierr;
  
    Range added_tets;
    for(_IT_CUBITMESHSETS_BY_BCDATA_TYPE_FOR_LOOP_(mField,BLOCKSET|BODYFORCESSET,it)) {
      int id = it->get_msId();
      EntityHandle meshset = it->get_meshset();
      rval = mField.get_moab().get_entities_by_type(meshset,MBTET,setOfBlocks[id].tEts,true); CHKERR_PETSC(rval);
      added_tets.merge(setOfBlocks[id].tEts);
      Block_BodyForces mydata;
      ierr = it->get_attribute_data_structure(mydata); CHKERRQ(ierr);
      setOfBlocks[id].rho0 = mydata.data.density;
      setOfBlocks[id].a0.resize(3);
      setOfBlocks[id].a0[0] = mydata.data.acceleration_x;
      setOfBlocks[id].a0[1] = mydata.data.acceleration_y;
      setOfBlocks[id].a0[2] = mydata.data.acceleration_z;
      //cerr << setOfBlocks[id].tEts << endl;
    }

    for(_IT_CUBITMESHSETS_BY_BCDATA_TYPE_FOR_LOOP_(mField,BLOCKSET|MAT_ELASTICSET,it)) {
      Mat_Elastic mydata;
      ierr = it->get_attribute_data_structure(mydata); CHKERRQ(ierr);
      if(mydata.data.User1 == 0) continue;
      Range tets;
      EntityHandle meshset = it->get_meshset();
      rval = mField.get_moab().get_entities_by_type(meshset,MBTET,tets,true); CHKERR_PETSC(rval);
      tets = subtract(tets,added_tets);
      if(tets.empty()) continue;
      int id = it->get_msId();
      setOfBlocks[-id].tEts = tets;
      setOfBlocks[-id].rho0 = mydata.data.User1;
      setOfBlocks[-id].a0.resize(3);
      setOfBlocks[-id].a0[0] = mydata.data.User2;
      setOfBlocks[-id].a0[1] = mydata.data.User3;
      setOfBlocks[-id].a0[2] = mydata.data.User4;
      //cerr << setOfBlocks[id].tEts << endl;
    }


    PetscFunctionReturn(0);
  }

  PetscErrorCode addConvectiveMassElement(string element_name,
    string velocity_field_name,
    string spatial_position_field_name,
    string material_position_field_name = "MESH_NODE_POSITIONS",
    bool ale = false,BitRefLevel bit = BitRefLevel()) {
    PetscFunctionBegin;

    PetscErrorCode ierr;
    //ErrorCode rval;

    ierr = mField.add_finite_element(element_name,MF_ZERO); CHKERRQ(ierr);
    ierr = mField.modify_finite_element_add_field_col(element_name,velocity_field_name); CHKERRQ(ierr);
    ierr = mField.modify_finite_element_add_field_data(element_name,velocity_field_name); CHKERRQ(ierr);
    ierr = mField.modify_finite_element_add_field_row(element_name,spatial_position_field_name); CHKERRQ(ierr);
    ierr = mField.modify_finite_element_add_field_col(element_name,spatial_position_field_name); CHKERRQ(ierr);
    ierr = mField.modify_finite_element_add_field_data(element_name,spatial_position_field_name); CHKERRQ(ierr);
    if(mField.check_field(material_position_field_name)) {
      if(ale) {
	ierr = mField.modify_finite_element_add_field_col(element_name,material_position_field_name); CHKERRQ(ierr);
	ierr = mField.modify_finite_element_add_field_data(element_name,"DOT_"+material_position_field_name); CHKERRQ(ierr);
      }
      ierr = mField.modify_finite_element_add_field_data(element_name,material_position_field_name); CHKERRQ(ierr);
    }
    ierr = mField.modify_finite_element_add_field_data(element_name,"DOT_"+velocity_field_name); CHKERRQ(ierr);
    ierr = mField.modify_finite_element_add_field_data(element_name,"DOT_"+spatial_position_field_name); CHKERRQ(ierr);

    Range tets;
    if(bit.any()) {
      ierr = mField.get_entities_by_type_and_ref_level(bit,BitRefLevel().set(),MBTET,tets); CHKERRQ(ierr);
    }

    map<int,BlockData>::iterator sit = setOfBlocks.begin();
    for(;sit!=setOfBlocks.end();sit++) {
      Range add_tets = sit->second.tEts;
      if(!tets.empty()) {
	add_tets = intersect(add_tets,tets);
      }
      ierr = mField.add_ents_to_finite_element_by_TETs(add_tets,element_name); CHKERRQ(ierr);
    }

    PetscFunctionReturn(0);
  }

  PetscErrorCode addVelocityElement(string element_name,
    string velocity_field_name,
    string spatial_position_field_name,
    string material_position_field_name = "MESH_NODE_POSITIONS",
    bool ale = false,BitRefLevel bit = BitRefLevel()) {
    PetscFunctionBegin;

    PetscErrorCode ierr;
    //ErrorCode rval;

    ierr = mField.add_finite_element(element_name,MF_ZERO); CHKERRQ(ierr);
    ierr = mField.modify_finite_element_add_field_row(element_name,velocity_field_name); CHKERRQ(ierr);
    ierr = mField.modify_finite_element_add_field_col(element_name,velocity_field_name); CHKERRQ(ierr);
    ierr = mField.modify_finite_element_add_field_data(element_name,velocity_field_name); CHKERRQ(ierr);
    ierr = mField.modify_finite_element_add_field_col(element_name,spatial_position_field_name); CHKERRQ(ierr);
    ierr = mField.modify_finite_element_add_field_data(element_name,spatial_position_field_name); CHKERRQ(ierr);
    if(mField.check_field(material_position_field_name)) {
      if(ale) {
	ierr = mField.modify_finite_element_add_field_col(element_name,material_position_field_name); CHKERRQ(ierr);
	ierr = mField.modify_finite_element_add_field_data(element_name,"DOT_"+material_position_field_name); CHKERRQ(ierr);
      }
      ierr = mField.modify_finite_element_add_field_data(element_name,material_position_field_name); CHKERRQ(ierr);
    }
    ierr = mField.modify_finite_element_add_field_data(element_name,"DOT_"+velocity_field_name); CHKERRQ(ierr);
    ierr = mField.modify_finite_element_add_field_data(element_name,"DOT_"+spatial_position_field_name); CHKERRQ(ierr);

    Range tets;
    if(bit.any()) {
      ierr = mField.get_entities_by_type_and_ref_level(bit,BitRefLevel().set(),MBTET,tets); CHKERRQ(ierr);
    }

    map<int,BlockData>::iterator sit = setOfBlocks.begin();
    for(;sit!=setOfBlocks.end();sit++) {
      Range add_tets = sit->second.tEts;
      if(!tets.empty()) {
	add_tets = intersect(add_tets,tets);
      }
      ierr = mField.add_ents_to_finite_element_by_TETs(add_tets,element_name); CHKERRQ(ierr);
    }

    PetscFunctionReturn(0);
  }


  PetscErrorCode addEshelbyDynamicMaterialMomentum(string element_name,
    string velocity_field_name,
    string spatial_position_field_name,
    string material_position_field_name = "MESH_NODE_POSITIONS",
    bool ale = false,BitRefLevel bit = BitRefLevel()) {
    PetscFunctionBegin;

    PetscErrorCode ierr;
    //ErrorCode rval;

    ierr = mField.add_finite_element(element_name,MF_ZERO); CHKERRQ(ierr);
    ierr = mField.modify_finite_element_add_field_col(element_name,velocity_field_name); CHKERRQ(ierr);
    ierr = mField.modify_finite_element_add_field_data(element_name,velocity_field_name); CHKERRQ(ierr);
    ierr = mField.modify_finite_element_add_field_col(element_name,spatial_position_field_name); CHKERRQ(ierr);
    ierr = mField.modify_finite_element_add_field_data(element_name,spatial_position_field_name); CHKERRQ(ierr);
    if(mField.check_field(material_position_field_name)) {
      if(ale) {
	ierr = mField.modify_finite_element_add_field_row(element_name,material_position_field_name); CHKERRQ(ierr);
	ierr = mField.modify_finite_element_add_field_col(element_name,material_position_field_name); CHKERRQ(ierr);
	ierr = mField.modify_finite_element_add_field_data(element_name,"DOT_"+material_position_field_name); CHKERRQ(ierr);
      }
      ierr = mField.modify_finite_element_add_field_data(element_name,material_position_field_name); CHKERRQ(ierr);
    }
    ierr = mField.modify_finite_element_add_field_data(element_name,"DOT_"+velocity_field_name); CHKERRQ(ierr);
    ierr = mField.modify_finite_element_add_field_data(element_name,"DOT_"+spatial_position_field_name); CHKERRQ(ierr);

    Range tets;
    if(bit.any()) {
      ierr = mField.get_entities_by_type_and_ref_level(bit,BitRefLevel().set(),MBTET,tets); CHKERRQ(ierr);
    }

    map<int,BlockData>::iterator sit = setOfBlocks.begin();
    for(;sit!=setOfBlocks.end();sit++) {
      Range add_tets = sit->second.tEts;
      if(!tets.empty()) {
	add_tets = intersect(add_tets,tets);
      }
      ierr = mField.add_ents_to_finite_element_by_TETs(add_tets,element_name); CHKERRQ(ierr);
    }

    PetscFunctionReturn(0);
  }

  PetscErrorCode setConvectiveMassOperators(
    string velocity_field_name,
    string spatial_position_field_name,
    string material_position_field_name = "MESH_NODE_POSITIONS",
    bool ale = false) {
    PetscFunctionBegin;

    commonData.spatialPositions = spatial_position_field_name;
    commonData.meshPositions = material_position_field_name;
    commonData.spatialVelocities = velocity_field_name;

    //Rhs
    feMassRhs.get_op_to_do_Rhs().push_back(new OpGetCommonDataAtGaussPts(velocity_field_name,commonData));
    feMassRhs.get_op_to_do_Rhs().push_back(new OpGetCommonDataAtGaussPts(spatial_position_field_name,commonData));
    feMassRhs.get_op_to_do_Rhs().push_back(new OpGetCommonDataAtGaussPts("DOT_"+velocity_field_name,commonData));
    feMassRhs.get_op_to_do_Rhs().push_back(new OpGetCommonDataAtGaussPts("DOT_"+spatial_position_field_name,commonData));
    if(mField.check_field(material_position_field_name)) {
      feMassRhs.get_op_to_do_Rhs().push_back(new OpGetCommonDataAtGaussPts(material_position_field_name,commonData));
      if(ale) {
	feMassRhs.get_op_to_do_Rhs().push_back(new OpGetCommonDataAtGaussPts("DOT_"+material_position_field_name,commonData));
      }
    }
    map<int,BlockData>::iterator sit = setOfBlocks.begin();
    for(;sit!=setOfBlocks.end();sit++) {
      feMassRhs.get_op_to_do_Rhs().push_back(new OpMassJacobian(spatial_position_field_name,sit->second,commonData,tAg,false));
      feMassRhs.get_op_to_do_Rhs().push_back(new OpMassRhs(spatial_position_field_name,sit->second,commonData));
    }

    //Lhs
    feMassLhs.get_op_to_do_Rhs().push_back(new OpGetCommonDataAtGaussPts(velocity_field_name,commonData));
    feMassLhs.get_op_to_do_Rhs().push_back(new OpGetCommonDataAtGaussPts(spatial_position_field_name,commonData));
    feMassLhs.get_op_to_do_Rhs().push_back(new OpGetCommonDataAtGaussPts("DOT_"+velocity_field_name,commonData));
    feMassLhs.get_op_to_do_Rhs().push_back(new OpGetCommonDataAtGaussPts("DOT_"+spatial_position_field_name,commonData));
    if(mField.check_field(material_position_field_name)) {
      feMassLhs.get_op_to_do_Rhs().push_back(new OpGetCommonDataAtGaussPts(material_position_field_name,commonData));
      if(ale) {
	feMassLhs.get_op_to_do_Rhs().push_back(new OpGetCommonDataAtGaussPts("DOT_"+material_position_field_name,commonData));
      }
    }
    sit = setOfBlocks.begin();
    for(;sit!=setOfBlocks.end();sit++) {
      feMassLhs.get_op_to_do_Rhs().push_back(new OpMassJacobian(spatial_position_field_name,sit->second,commonData,tAg));
      feMassLhs.get_op_to_do_Lhs().push_back(new OpMassLhs_dM_dv(spatial_position_field_name,velocity_field_name,sit->second,commonData));
      feMassLhs.get_op_to_do_Lhs().push_back(new OpMassLhs_dM_dx(spatial_position_field_name,spatial_position_field_name,sit->second,commonData));
      if(mField.check_field(material_position_field_name)) {
	if(ale) {
	  feMassLhs.get_op_to_do_Lhs().push_back(new OpMassLhs_dM_dX(spatial_position_field_name,material_position_field_name,sit->second,commonData));
	}
      }
    }

    PetscFunctionReturn(0);
  }

  PetscErrorCode setVelocityOperators(
    string velocity_field_name,
    string spatial_position_field_name,
    string material_position_field_name = "MESH_NODE_POSITIONS",
    bool ale = false) {
    PetscFunctionBegin;

    commonData.spatialPositions = spatial_position_field_name;
    commonData.meshPositions = material_position_field_name;
    commonData.spatialVelocities = velocity_field_name;

    //Rhs
    feVelRhs.get_op_to_do_Rhs().push_back(new OpGetCommonDataAtGaussPts(velocity_field_name,commonData));
    feVelRhs.get_op_to_do_Rhs().push_back(new OpGetCommonDataAtGaussPts(spatial_position_field_name,commonData));
    feVelRhs.get_op_to_do_Rhs().push_back(new OpGetCommonDataAtGaussPts("DOT_"+velocity_field_name,commonData));
    feVelRhs.get_op_to_do_Rhs().push_back(new OpGetCommonDataAtGaussPts("DOT_"+spatial_position_field_name,commonData));
    if(mField.check_field(material_position_field_name)) {
      feVelRhs.get_op_to_do_Rhs().push_back(new OpGetCommonDataAtGaussPts(material_position_field_name,commonData));
      if(ale) {
	feVelRhs.get_op_to_do_Rhs().push_back(new OpGetCommonDataAtGaussPts("DOT_"+material_position_field_name,commonData));
      }
    }
    map<int,BlockData>::iterator sit = setOfBlocks.begin();
    for(;sit!=setOfBlocks.end();sit++) {
      feVelRhs.get_op_to_do_Rhs().push_back(new OpVelocityJacobian(velocity_field_name,sit->second,commonData,tAg,false));
      feVelRhs.get_op_to_do_Rhs().push_back(new OpVelocityRhs(velocity_field_name,sit->second,commonData));
    }

    //Lhs
    feVelLhs.get_op_to_do_Rhs().push_back(new OpGetCommonDataAtGaussPts(velocity_field_name,commonData));
    feVelLhs.get_op_to_do_Rhs().push_back(new OpGetCommonDataAtGaussPts(spatial_position_field_name,commonData));
    feVelLhs.get_op_to_do_Rhs().push_back(new OpGetCommonDataAtGaussPts("DOT_"+velocity_field_name,commonData));
    feVelLhs.get_op_to_do_Rhs().push_back(new OpGetCommonDataAtGaussPts("DOT_"+spatial_position_field_name,commonData));
    if(mField.check_field(material_position_field_name)) {
      feVelLhs.get_op_to_do_Rhs().push_back(new OpGetCommonDataAtGaussPts(material_position_field_name,commonData));
      if(ale) {
	feVelLhs.get_op_to_do_Rhs().push_back(new OpGetCommonDataAtGaussPts("DOT_"+material_position_field_name,commonData));
      }
    }
    sit = setOfBlocks.begin();
    for(;sit!=setOfBlocks.end();sit++) {
      feVelLhs.get_op_to_do_Rhs().push_back(new OpVelocityJacobian(velocity_field_name,sit->second,commonData,tAg));
      feVelLhs.get_op_to_do_Lhs().push_back(new OpVelocityLhs_dV_dv(velocity_field_name,velocity_field_name,sit->second,commonData));
      feVelLhs.get_op_to_do_Lhs().push_back(new OpVelocityLhs_dV_dx(velocity_field_name,spatial_position_field_name,sit->second,commonData));
      if(mField.check_field(material_position_field_name)) {
	if(ale) {
	  feVelLhs.get_op_to_do_Lhs().push_back(new OpVelocityLhs_dV_dX(velocity_field_name,material_position_field_name,sit->second,commonData));
	}
      }
    }

    PetscFunctionReturn(0);
  }


  PetscErrorCode setKinematicEshelbyOperators(
    string velocity_field_name,
    string spatial_position_field_name,
    string material_position_field_name = "MESH_NODE_POSITIONS",bool ale = false) {
    PetscFunctionBegin;

    commonData.spatialPositions = spatial_position_field_name;
    commonData.meshPositions = material_position_field_name;
    commonData.spatialVelocities = velocity_field_name;

    //Rhs
    feTRhs.get_op_to_do_Rhs().push_back(new OpGetCommonDataAtGaussPts(velocity_field_name,commonData));
    feTRhs.get_op_to_do_Rhs().push_back(new OpGetCommonDataAtGaussPts(spatial_position_field_name,commonData));
    if(mField.check_field(material_position_field_name)) {
      feTRhs.get_op_to_do_Rhs().push_back(new OpGetCommonDataAtGaussPts(material_position_field_name,commonData));
    }
    map<int,BlockData>::iterator sit = setOfBlocks.begin();
    for(;sit!=setOfBlocks.end();sit++) {
      feTRhs.get_op_to_do_Rhs().push_back(new OpEshelbyDynamicMaterialMomentumJacobian(material_position_field_name,sit->second,commonData,tAg,false));
      feTRhs.get_op_to_do_Rhs().push_back(new OpEshelbyDynamicMaterialMomentumRhs(material_position_field_name,sit->second,commonData));
    }

    //Lhs
    feTLhs.get_op_to_do_Rhs().push_back(new OpGetCommonDataAtGaussPts(velocity_field_name,commonData));
    feTLhs.get_op_to_do_Rhs().push_back(new OpGetCommonDataAtGaussPts(spatial_position_field_name,commonData));
    if(mField.check_field(material_position_field_name)) {
      feTLhs.get_op_to_do_Rhs().push_back(new OpGetCommonDataAtGaussPts(material_position_field_name,commonData));
    }
    sit = setOfBlocks.begin();
    for(;sit!=setOfBlocks.end();sit++) {
      feTLhs.get_op_to_do_Rhs().push_back(new OpEshelbyDynamicMaterialMomentumJacobian(material_position_field_name,sit->second,commonData,tAg));
      feTLhs.get_op_to_do_Lhs().push_back(new OpEshelbyDynamicMaterialMomentumRhs_dv(material_position_field_name,velocity_field_name,sit->second,commonData));
      feTLhs.get_op_to_do_Lhs().push_back(new OpEshelbyDynamicMaterialMomentumRhs_dx(material_position_field_name,spatial_position_field_name,sit->second,commonData));
      if(mField.check_field(material_position_field_name)) {
	if(ale) {
	  feTLhs.get_op_to_do_Lhs().push_back(new OpEshelbyDynamicMaterialMomentumRhs_dX(material_position_field_name,material_position_field_name,sit->second,commonData));
	}
      }
    }

    PetscFunctionReturn(0);
  }

};


#endif //__CONVECTIVE_MASS_ELEMENT_HPP

/***************************************************************************//**
 * \defgroup convective_mass_elem (Convective) Mass Element
 * \ingroup mofem_forces_and_sources 
 ******************************************************************************/


