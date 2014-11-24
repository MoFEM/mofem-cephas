/** \file ThermalElement.hpp 
 * \brief Operators and data structures for thermal analys
 *
 * Implementation of thermal element for unsteady and steady case.
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

  FieldInterface &mField;
  short int tAg;

  ConvectiveMassElement(
    FieldInterface &m_field,short int tag):
    feMassRhs(m_field),feMassLhs(m_field),
    feVelRhs(m_field),feVelLhs(m_field),
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
    PetscErrorCode iNvert(ublas::matrix<TYPE> a,ublas::matrix<TYPE> &inv_a) {
      PetscFunctionBegin;
      PetscErrorCode ierr;
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
      TYPE det;
      ierr = dEterminatnt(a,det); CHKERRQ(ierr);
      inv_a /= det;
      PetscFunctionReturn(0);
    }
  
    template<typename TYPE> 
    PetscErrorCode calculateMomentumRate(
      double rho0,ublas::vector<double>& a0,
      ublas::vector<TYPE>& a,ublas::matrix<TYPE>& grad_v,
      ublas::vector<TYPE>& dot_W,ublas::matrix<TYPE>& H,
      ublas::matrix<TYPE>& h,
      ublas::vector<TYPE>& dp_dt) {
      PetscFunctionBegin;
  
      PetscErrorCode ierr;
  
      //calulate gradient of deformation
      ublas::matrix<TYPE> invH(3,3);
      ierr = iNvert(H,invH); CHKERRQ(ierr);
      ublas::matrix<TYPE> F(3,3);
      noalias(F) = prod(h,invH);
      TYPE detF;
      ierr = dEterminatnt(F,detF); CHKERRQ(ierr);
  
      //calulate current density
      TYPE rho = rho0*detF*rho0;
  
      //momentum rate
      noalias(dp_dt) = rho*(a0 + a + prod(grad_v,dot_W));
  
      PetscFunctionReturn(0);
    }
  
    template<typename TYPE> 
    PetscErrorCode calculateFuncUnderIntegral(
      double rho0,ublas::vector<double>& a0,
      ublas::vector<TYPE>& a,
      ublas::matrix<TYPE>& grad_v,
      ublas::vector<TYPE>& dot_W,
      ublas::matrix<TYPE>& H,
      ublas::matrix<TYPE>& h,
      ublas::vector<TYPE>& f) {
      PetscFunctionBegin;
  
      PetscErrorCode ierr;
      ublas::vector<TYPE> dp_dt;
      dp_dt.resize(3);
      ierr = calculateMomentumRate(
        rho0,a0,a,grad_v,dot_W,H,h,dp_dt); CHKERRQ(ierr);
      TYPE detH;
      ierr = dEterminatnt(H,detH); CHKERRQ(ierr);
      noalias(f) = dp_dt*detH;

      PetscFunctionReturn(0);
    }

    template<typename TYPE> 
    PetscErrorCode calulateVelocity(
      ublas::vector<TYPE>& dot_w,
      ublas::vector<TYPE>& dot_W,
      ublas::matrix<TYPE>& h,
      ublas::matrix<TYPE>& H,
      ublas::vector<TYPE>& dot_u) {
      PetscFunctionBegin;
      PetscErrorCode ierr;

      ublas::matrix<TYPE> invH;
      ierr = iNvert(H,invH); CHKERRQ(ierr);
      ublas::matrix<TYPE> F;
      F = prod(h,invH);

      dot_u = dot_w + prod(F,dot_W);

      PetscFunctionReturn(0);
    }

  };

  struct OpMassRhs: public TetElementForcesAndSourcesCore::UserDataOperator,CommonFunctions {

    BlockData &dAta;
    CommonData &commonData;

    OpMassRhs(const string field_name,BlockData &data,CommonData &common_data):
      TetElementForcesAndSourcesCore::UserDataOperator(field_name),
      dAta(data),commonData(common_data) { }

    ublas::vector<double> a;
    ublas::vector<double> nf;
    ublas::vector<double> dot_W;
    ublas::matrix<double> H;
    ublas::matrix<double> h;
    ublas::matrix<double> g;
    ublas::vector<double> f;
 
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

	a.resize(3);
	dot_W.resize(3);
	H.resize(3,3);
	h.resize(3,3);
	g.resize(3,3);
	f.resize(3);

	dot_W.clear();
	H.clear();
	for(int dd = 0;dd<3;dd++) {
	  H(dd,dd) = 1;
	}

	nf.resize(nb_dofs);
	nf.clear();

	for(unsigned int gg = 0;gg<row_data.getN().size1();gg++) {

	  noalias(a) = commonData.dataAtGaussPts["DOT_"+commonData.spatialVelocities][gg];
	  noalias(g) = commonData.gradAtGaussPts[commonData.spatialVelocities][gg];
	  noalias(h) = commonData.gradAtGaussPts[commonData.spatialPositions][gg];

	  if(commonData.dataAtGaussPts["DOT_"+commonData.meshPositions].size()>0) {
	    noalias(dot_W) = commonData.dataAtGaussPts["DOT_"+commonData.meshPositions][gg];
	  } 
	  if(commonData.gradAtGaussPts[commonData.meshPositions].size()>0) {
	    noalias(H) = commonData.gradAtGaussPts[commonData.meshPositions][gg];
	  } 

	  double rho0 = dAta.rho0;
	  ublas::vector<double>& a0 = dAta.a0;
	  ierr = calculateFuncUnderIntegral(
	    rho0,a0,a,g,dot_W,H,h,f); CHKERRQ(ierr);
	  double val = getVolume()*getGaussPts()(3,gg);
	  f *= val;

	  //cerr << f << endl;

	  for(unsigned int dd = 0;dd<nb_dofs/3;dd++) {
	    for(int rr = 0;rr<3;rr++) {
	      nf[3*dd+rr] += row_data.getN()(gg,dd)*f[rr];
	    }
	  }

	}

	if(nb_dofs > 3*row_data.getN().size2()) {
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

  struct OpMassLhs_dM_dX: public TetElementForcesAndSourcesCore::UserDataOperator,CommonFunctions {

    BlockData &dAta;
    CommonData &commonData;
    int tAg;

    OpMassLhs_dM_dX(const string field_name,const string col_field,BlockData &data,CommonData &common_data,int tag):
      TetElementForcesAndSourcesCore::UserDataOperator(field_name,col_field),
      dAta(data),commonData(common_data),tAg(tag) { symm = false; }

    ublas::vector<adouble> a;
    ublas::matrix<adouble> g;
    ublas::vector<adouble> dot_W;
    ublas::matrix<adouble> H;
    ublas::matrix<adouble> h;

    ublas::vector<double> dX;
    ublas::vector<adouble> a_dX;
    ublas::vector<adouble> a_f;
    double *active_ptr;
    int nb_active_vars;

    virtual PetscErrorCode setActive(DataForcesAndSurcesCore::EntData &col_data,int gg) {
      PetscFunctionBegin;

      //cerr << "OpMassLhs_dM_dX" << endl;
      
      //active
      int nb_dofs = col_data.getIndices().size();
      dX.resize(nb_dofs);
      dX.clear();
      a_dX.resize(nb_dofs);
      for(unsigned int nn = 0;nn<nb_dofs;nn++) {
	a_dX[nn] <<= dX[nn];
      }
      ublas::vector<double> N = col_data.getN(gg,nb_dofs/3);
      ublas::matrix<double> diffN = col_data.getDiffN(gg,nb_dofs/3);
      for(unsigned int nn1 = 0;nn1<3;nn1++) {
	for(unsigned int dd = 0;dd<nb_dofs/3;dd++) {
	  dot_W[nn1] += N[dd]*a_dX[3*dd+nn1]*getFEMethod()->ts_a;
	}
	for(unsigned int nn2 = 0;nn2<3;nn2++) {
	  for(unsigned int dd = 0;dd<nb_dofs/3;dd++) {
	    H(nn1,nn2) += diffN(dd,nn2)*a_dX[3*dd+nn1];
	  }
	}
      }
      active_ptr = &*dX.data().begin();
      nb_active_vars = dX.size();

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
      if(row_data.getIndices().size()==0) PetscFunctionReturn(0);
      if(col_data.getIndices().size()==0) PetscFunctionReturn(0);

      int nb_row = row_data.getIndices().size();
      int nb_col = col_data.getIndices().size();

      try {

	vector<double*> jac_row_ptr;
	ublas::vector<double> f;
	ublas::matrix<double> jac;
	ublas::matrix<double> k;
	a_f.resize(3);

	a.resize(3);
	h.resize(3,3);
	g.resize(3,3);

	dot_W.resize(3);
	dot_W.clear();
	H.resize(3,3);
	H.clear();
	for(int dd = 0;dd<3;dd++) {
	  H(dd,dd) = 1;
	}

	k.resize(nb_row,nb_col);
	k.clear();

	for(unsigned int gg = 0;gg<row_data.getN().size1();gg++) {

	  //set passive variables
	  noalias(a) = commonData.dataAtGaussPts["DOT_"+commonData.spatialVelocities][gg];
	  noalias(g) = commonData.gradAtGaussPts[commonData.spatialVelocities][gg];
	  noalias(h) = commonData.gradAtGaussPts[commonData.spatialPositions][gg];

	  if(commonData.dataAtGaussPts["DOT_"+commonData.meshPositions].size()>0) {
	    noalias(dot_W) = commonData.dataAtGaussPts["DOT_"+commonData.meshPositions][gg];
	  } 
	  if(commonData.gradAtGaussPts[commonData.meshPositions].size()>0) {
	    noalias(H) = commonData.gradAtGaussPts[commonData.meshPositions][gg];
	  } 

	  trace_on(tAg);

	  //set active variables
	  ierr = setActive(col_data,gg); CHKERRQ(ierr);

	  double rho0 = dAta.rho0;
	  ublas::vector<double>& a0 = dAta.a0;
	  ierr = calculateFuncUnderIntegral(
	    rho0,a0,a,g,dot_W,H,h,a_f); CHKERRQ(ierr);
	  double val = getVolume()*getGaussPts()(3,gg);
	  a_f *= val;

	  //dependant
	  f.resize(3);
	  for(int rr = 0;rr<3;rr++) {
	    a_f[rr] >>= f[rr];
	  }

	  trace_off();

	  if(gg == 0) {
	    jac_row_ptr.resize(3);
	    jac.resize(3,nb_active_vars);
	    for(int nn1 = 0;nn1<3;nn1++) {
	      jac_row_ptr[nn1] = &jac(nn1,0);     
	    }
	    if(nb_active_vars!=nb_col) {
	      SETERRQ2(PETSC_COMM_SELF,1,"data inconsistency %d!=%d",nb_active_vars,col_data.getIndices().size());
	    }
	  }

	  int r;
	  jac.clear();
	  r = jacobian(
	    tAg,3,nb_active_vars,
	    active_ptr,&jac_row_ptr[0]);
	  //cerr << "jac " << jac << endl;
	  //cerr << row_data.getIndices() << endl;
	  //cerr << col_data.getIndices() << endl;

	  { //integrate element stiffnes matrix
	    for(unsigned int dd1 = 0;dd1<nb_row/3;dd1++) {
	      for(int rr1 = 0;rr1<3;rr1++) {
		for(unsigned int dd2 = 0;dd2<nb_col/3;dd2++) {
		  for(int rr2 = 0;rr2<3;rr2++) {
		    k(3*dd1+rr1,3*dd2+rr2) += row_data.getN()(gg,dd1)*jac(rr1,3*dd2+rr2);
		  }
		}
	      }
	    }
	  }
	  //cerr << k << endl;

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

  struct OpMassLhs_dM_dx: public OpMassLhs_dM_dX  {

    OpMassLhs_dM_dx(const string field_name,const string col_field,BlockData &data,CommonData &common_data,int tag):
      OpMassLhs_dM_dX(field_name,col_field,data,common_data,tag) {}

    ublas::vector<double> dx;
    ublas::vector<adouble> a_dx;

    PetscErrorCode setActive(DataForcesAndSurcesCore::EntData &col_data,int gg) {
      PetscFunctionBegin;

      //cerr << "OpMassLhs_dM_dx" << endl;
      
      //active
      int nb_dofs = col_data.getIndices().size();
      dx.resize(nb_dofs);
      dx.clear();
      a_dx.resize(nb_dofs);
      for(unsigned int nn = 0;nn<nb_dofs;nn++) {
	a_dx[nn] <<= dx[nn];
      }
      ublas::matrix<double> diffN = col_data.getDiffN(gg,nb_dofs/3);
      for(unsigned int nn1 = 0;nn1<3;nn1++) {
	for(unsigned int nn2 = 0;nn2<3;nn2++) {
	  for(unsigned int dd = 0;dd<nb_dofs/3;dd++) {
	    h(nn1,nn2) += diffN(dd,nn2)*a_dx[3*dd+nn1];
	  }
	}
      }
      active_ptr = &*dx.data().begin();
      nb_active_vars = nb_dofs;

      PetscFunctionReturn(0);
    }

  };

  struct OpMassLhs_dM_dv: public OpMassLhs_dM_dX  {

    OpMassLhs_dM_dv(const string field_name,const string col_field,BlockData &data,CommonData &common_data,int tag):
      OpMassLhs_dM_dX(field_name,col_field,data,common_data,tag) {}

    ublas::vector<double> dv;
    ublas::vector<adouble> a_dv;

    PetscErrorCode setActive(DataForcesAndSurcesCore::EntData &col_data,int gg) {
      PetscFunctionBegin;
    
      //cerr << "OpMassLhs_dM_dv" << endl;

      try {

      //active
      int nb_dofs = col_data.getIndices().size();
      dv.resize(nb_dofs);
      dv.clear();
      a_dv.resize(nb_dofs);
      for(unsigned int nn = 0;nn<nb_dofs;nn++) {
	a_dv[nn] <<= dv[nn];
      }

      ublas::vector<double> N = col_data.getN(gg,nb_dofs/3);
      ublas::matrix<double> diffN = col_data.getDiffN(gg,nb_dofs/3);
      for(unsigned int nn1 = 0;nn1<3;nn1++) {
	for(unsigned int dd = 0;dd<nb_dofs/3;dd++) {
	  a[nn1] += N[dd]*a_dv[3*dd+nn1]*getFEMethod()->ts_a;
	}
	for(unsigned int nn2 = 0;nn2<3;nn2++) {
	  for(unsigned int dd = 0;dd<nb_dofs/3;dd++) {
	    g(nn1,nn2) += diffN(dd,nn2)*a_dv[3*dd+nn1];
	  }
	}
      }
      active_ptr = &*dv.data().begin();
      nb_active_vars = nb_dofs;

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

    ublas::vector<double> dot_W;
    ublas::matrix<double> H;
    ublas::vector<double> dot_u;
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

	dot_W.resize(3);
	dot_W.clear();
	H.resize(3,3);
	H.clear();
	dot_u.resize(3);
	for(int dd = 0;dd<3;dd++) {
	  H(dd,dd) = 1;
	}

	nf.resize(nb_dofs,0);
	nf.clear();
  
	for(unsigned int gg = 0;gg<row_data.getN().size1();gg++) {

	  ublas::vector<double>& v = commonData.dataAtGaussPts[commonData.spatialVelocities][gg];
	  ublas::vector<double>& dot_w = commonData.dataAtGaussPts["DOT_"+commonData.spatialPositions][gg];
	  ublas::matrix<double>& h = commonData.gradAtGaussPts[commonData.spatialPositions][gg];
	  if(commonData.dataAtGaussPts["DOT_"+commonData.meshPositions].size()>0) {
	   noalias(dot_W) = commonData.dataAtGaussPts["DOT_"+commonData.meshPositions][gg];
	  }
	  if(commonData.gradAtGaussPts[commonData.meshPositions].size()>0) {
	    noalias(H) = commonData.gradAtGaussPts[commonData.meshPositions][gg];
	  }

	  ierr = calulateVelocity(dot_w,dot_W,h,H,dot_u); CHKERRQ(ierr);
	  double detH = 1;
	  if(commonData.gradAtGaussPts[commonData.meshPositions].size()>0) {
	    ierr = dEterminatnt(H,detH); CHKERRQ(ierr);
	  }
	  ublas::vector<double> res = (v - dot_u)*detH;

	  double val = getVolume()*getGaussPts()(3,gg);
	  res *= val;

	  for(unsigned int dd = 0;dd<nb_dofs/3;dd++) {
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

  struct OpVelocityLhs_Jacobian: public TetElementForcesAndSourcesCore::UserDataOperator,CommonFunctions {

    BlockData &dAta;
    CommonData &commonData;
    int tAg;

    OpVelocityLhs_Jacobian(const string vel_field,const string field_name,BlockData &data,CommonData &common_data,int tag):
      TetElementForcesAndSourcesCore::UserDataOperator(vel_field,field_name),
      dAta(data),commonData(common_data),tAg(tag) { symm = false;  }


  };


  struct OpVelocityLhs_dV_dX: public TetElementForcesAndSourcesCore::UserDataOperator,CommonFunctions {

    BlockData &dAta;
    CommonData &commonData;
    int tAg;

    OpVelocityLhs_dV_dX(const string vel_field,const string field_name,BlockData &data,CommonData &common_data,int tag):
      TetElementForcesAndSourcesCore::UserDataOperator(vel_field,field_name),
      dAta(data),commonData(common_data),tAg(tag) { symm = false;  }

    ublas::vector<adouble> v;
    ublas::vector<adouble> dot_u;
    ublas::vector<adouble> dot_w;
    ublas::vector<adouble> dot_W;
    ublas::matrix<adouble> H;
    ublas::matrix<adouble> h;

    vector<double*> jac_row_ptr;
    ublas::matrix<double> jac;
    ublas::matrix<double> k;
    ublas::vector<double> res;
    ublas::vector<adouble> a_res;

    ublas::vector<double> dX;
    ublas::vector<adouble> a_dX;
    double *active_ptr;
    int nb_active_vars;

    virtual PetscErrorCode setActive(DataForcesAndSurcesCore::EntData &col_data,int gg) {
      PetscFunctionBegin;
      //cerr << "OpVelocityLhs_dV_dX" << endl;
      //active
      dX.resize(col_data.getIndices().size());
      dX.clear();
      a_dX.resize(dX.size());
      for(unsigned int nn = 0;nn<dX.size();nn++) {
	a_dX[nn] <<= dX[nn];
      }
      ublas::vector<double> N = col_data.getN(gg,dX.size()/3);
      ublas::matrix<double> diffN = trans(col_data.getDiffN(gg,a_dX.size()));
      for(unsigned int nn1 = 0;nn1<3;nn1++) {
	for(unsigned int dd = 0;dd<a_dX.size();dd++) {
	  dot_W[nn1] += N[dd]*a_dX[3*dd+nn1]*getFEMethod()->ts_a;
	}
	for(unsigned int nn2 = 0;nn2<3;nn2++) {
	  for(unsigned int dd = 0;dd<a_dX.size();dd++) {
	    H(nn1,nn2) += diffN(nn2,dd)*a_dX[3*dd+nn1];
	  }
	}
      }
      active_ptr = &*dX.data().begin();
      nb_active_vars = dX.size();
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

      //cerr << row_side << " " << col_side << " " << row_type << " " << col_type << " : " << nb_row << " " << nb_col << endl;

      try {

	v.resize(3);
	dot_w.resize(3);
	h.resize(3,3);
	dot_W.resize(3);
	dot_W.clear();
	H.resize(3,3);
	H.clear();
	for(int dd = 0;dd<3;dd++) {
	  H(dd,dd) = 1;
	}
	k.resize(nb_row,nb_col);
	k.clear();
	//cerr << k << endl;

	for(unsigned int gg = 0;gg<row_data.getN().size1();gg++) {
	  dot_u.resize(3);
	  a_res.resize(3);
	  //set active and passive variables
	  noalias(v) = commonData.dataAtGaussPts[commonData.spatialVelocities][gg];
	  noalias(dot_w) = commonData.dataAtGaussPts["DOT_"+commonData.spatialPositions][gg];
	  noalias(h) = commonData.gradAtGaussPts[commonData.spatialPositions][gg];
	  if(commonData.dataAtGaussPts["DOT_"+commonData.meshPositions].size()>0) {
	   noalias(dot_W) = commonData.dataAtGaussPts["DOT_"+commonData.meshPositions][gg];
	  }
	  if(commonData.gradAtGaussPts[commonData.meshPositions].size()>0) {
	    noalias(H) = commonData.gradAtGaussPts[commonData.meshPositions][gg];
	  }

	  //cerr << "K v: " << v << endl;
	  //cerr << "K dot_w: " << dot_w << endl;

	  trace_on(tAg);
	  ierr = setActive(col_data,gg); CHKERRQ(ierr);
	  ierr = calulateVelocity(dot_w,dot_W,h,H,dot_u); CHKERRQ(ierr);
	  //cerr << "dot_u " << dot_u << endl;
	  adouble detH;
	  ierr = dEterminatnt(H,detH); CHKERRQ(ierr);
	  noalias(a_res) = (v - dot_u)*detH;
	  //cerr << "a_res " << a_res << endl;
	  //dependant
	  res.resize(3);
	  for(int rr = 0;rr<3;rr++) {
	    a_res[rr] >>= res[rr];
	  }
	  trace_off();

	  //cerr << "res " << res << endl;

	  /*size_t tape_stats[11];
	  tapestats(tAg,tape_stats);
	  cerr << "the number of independents, i.e. calls to <<= " << tape_stats[0] << endl;
	  cerr << "the number of dependents, i.e. calls to >>= " << tape_stats[1] << endl;
	  cerr << "the maximal number of live active variables " << tape_stats[2] << endl;
	  cerr << "the size of value stack (number of overwrites) " << tape_stats[3] << endl;
	  cerr << "the buffer size (a multiple of eight) " << tape_stats[4] << endl;
	  cerr << "the total number of operations recorded " << tape_stats[5] << endl;*/

	  if(gg == 0) {
	    jac_row_ptr.resize(3);
	    jac.resize(3,nb_active_vars);
	    for(int nn1 = 0;nn1<3;nn1++) {
	      jac_row_ptr[nn1] = &jac(nn1,0);     
	    }
	    if(nb_active_vars!=nb_col) {
	      SETERRQ(PETSC_COMM_SELF,1,"data inconsistency");
	    }
	  }

	  int r;
	  jac.clear();
	  r = jacobian(
	    tAg,3,nb_active_vars,
	    active_ptr,&jac_row_ptr[0]);
	  double val = getVolume()*getGaussPts()(3,gg);
	  jac *= val;
	  //cerr << jac << endl;

	  { //integrate element stiffnes matrix
	    for(unsigned int dd1 = 0;dd1<nb_row/3;dd1++) {
	      for(int rr1 = 0;rr1<3;rr1++) {
		for(unsigned int dd2 = 0;dd2<nb_col/3;dd2++) {
		  for(int rr2 = 0;rr2<3;rr2++) {
		    k(3*dd1+rr1,3*dd2+rr2) += row_data.getN()(gg,dd1)*jac(rr1,3*dd2+rr2);
		  }
		}
	      }
	    }
	  }
	  //cerr << k << endl;
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

  struct OpVelocityLhs_dV_dx: public OpVelocityLhs_dV_dX {

    OpVelocityLhs_dV_dx(const string vel_field,const string field_name,BlockData &data,CommonData &common_data,int tag):
      OpVelocityLhs_dV_dX(vel_field,field_name,data,common_data,tag) {}

    ublas::vector<adouble> a_dx;
    ublas::vector<double> dx;

    PetscErrorCode setActive(DataForcesAndSurcesCore::EntData &col_data,int gg) {
      PetscFunctionBegin;

      //cerr << "OpVelocityLhs_dV_dx" << endl;

      //active
      int nb_dofs = col_data.getIndices().size();
      dx.resize(nb_dofs);
      dx.clear();
      a_dx.resize(nb_dofs);
      for(unsigned int nn = 0;nn<dx.size();nn++) {
	a_dx[nn] <<= dx[nn];
      }

      ublas::vector<double> N = col_data.getN(gg,nb_dofs/3);
      ublas::matrix<double> diffN = trans(col_data.getDiffN(gg,nb_dofs/3));
      for(unsigned int nn1 = 0;nn1<3;nn1++) {
	for(unsigned int dd = 0;dd<a_dx.size()/3;dd++) {
	  dot_w[nn1] += N[dd]*a_dx[3*dd+nn1]*getFEMethod()->ts_a;
	}
	for(unsigned int nn2 = 0;nn2<3;nn2++) {
	  for(unsigned int dd = 0;dd<a_dx.size()/3;dd++) {
	    h(nn1,nn2) += diffN(nn2,dd)*a_dx[3*dd+nn1];
	  }
	}
      }
      active_ptr = &*dx.data().begin();
      nb_active_vars = nb_dofs;

      PetscFunctionReturn(0);
    }

  };

  struct OpVelocityLhs_dV_dv: public OpVelocityLhs_dV_dX {

    OpVelocityLhs_dV_dv(const string vel_field,const string field_name,BlockData &data,CommonData &common_data,int tag):
      OpVelocityLhs_dV_dX(vel_field,field_name,data,common_data,tag) {}

    ublas::vector<adouble> a_dv;
    ublas::vector<double> dv;

    PetscErrorCode setActive(DataForcesAndSurcesCore::EntData &col_data,int gg) {
      PetscFunctionBegin;

      //cerr << "OpVelocityLhs_dV_dv" << endl;

      //active
      int nb_dofs = col_data.getIndices().size();
      dv.resize(nb_dofs);
      dv.clear();
      a_dv.resize(nb_dofs);
      for(unsigned int nn = 0;nn<nb_dofs;nn++) {
	a_dv[nn] <<= dv[nn];
      }

      ublas::vector<double> N = col_data.getN(gg,nb_dofs/3);
      for(unsigned int nn1 = 0;nn1<3;nn1++) {
	for(unsigned int dd = 0;dd<nb_dofs/3;dd++) {
	  v[nn1] += N[dd]*a_dv[3*dd+nn1];
	}
      }

      active_ptr = &*dv.data().begin();
      nb_active_vars = nb_dofs;

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
    string material_position_field_name = "MESH_NODE_POSITIONS",bool ale = false) {
    PetscFunctionBegin;

    PetscErrorCode ierr;
    //ErrorCode rval;

    ierr = mField.add_finite_element(element_name,MF_ZERO); CHKERRQ(ierr);
    ierr = mField.modify_finite_element_add_field_row(element_name,velocity_field_name); CHKERRQ(ierr);
    ierr = mField.modify_finite_element_add_field_col(element_name,velocity_field_name); CHKERRQ(ierr);
    ierr = mField.modify_finite_element_add_field_data(element_name,velocity_field_name); CHKERRQ(ierr);
    ierr = mField.modify_finite_element_add_field_row(element_name,spatial_position_field_name); CHKERRQ(ierr);
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

    map<int,BlockData>::iterator sit = setOfBlocks.begin();
    for(;sit!=setOfBlocks.end();sit++) {
      ierr = mField.add_ents_to_finite_element_by_TETs(sit->second.tEts,element_name); CHKERRQ(ierr);
    }

    PetscFunctionReturn(0);
  }

  PetscErrorCode addVelocityElement(string element_name,
    string velocity_field_name,
    string spatial_position_field_name,
    string material_position_field_name = "MESH_NODE_POSITIONS",bool ale = false) {
    PetscFunctionBegin;

    PetscErrorCode ierr;
    //ErrorCode rval;

    ierr = mField.add_finite_element(element_name,MF_ZERO); CHKERRQ(ierr);
    ierr = mField.modify_finite_element_add_field_row(element_name,velocity_field_name); CHKERRQ(ierr);
    ierr = mField.modify_finite_element_add_field_col(element_name,velocity_field_name); CHKERRQ(ierr);
    ierr = mField.modify_finite_element_add_field_data(element_name,velocity_field_name); CHKERRQ(ierr);
    ierr = mField.modify_finite_element_add_field_row(element_name,spatial_position_field_name); CHKERRQ(ierr);
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

    map<int,BlockData>::iterator sit = setOfBlocks.begin();
    for(;sit!=setOfBlocks.end();sit++) {
      mField.add_ents_to_finite_element_by_TETs(sit->second.tEts,element_name);
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
      feMassLhs.get_op_to_do_Lhs().push_back(new OpMassLhs_dM_dv(spatial_position_field_name,velocity_field_name,sit->second,commonData,tAg));
      feMassLhs.get_op_to_do_Lhs().push_back(new OpMassLhs_dM_dx(spatial_position_field_name,spatial_position_field_name,sit->second,commonData,tAg));
      if(mField.check_field(material_position_field_name)) {
	if(ale) {
	  feMassLhs.get_op_to_do_Lhs().push_back(new OpMassLhs_dM_dX(spatial_position_field_name,material_position_field_name,sit->second,commonData,tAg));
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
      feVelLhs.get_op_to_do_Lhs().push_back(new OpVelocityLhs_dV_dx(velocity_field_name,spatial_position_field_name,sit->second,commonData,tAg));
      feVelLhs.get_op_to_do_Lhs().push_back(new OpVelocityLhs_dV_dv(velocity_field_name,velocity_field_name,sit->second,commonData,tAg));
      if(mField.check_field(material_position_field_name)) {
	if(ale) {
	  feVelLhs.get_op_to_do_Lhs().push_back(new OpVelocityLhs_dV_dX(velocity_field_name,material_position_field_name,sit->second,commonData,tAg));
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



