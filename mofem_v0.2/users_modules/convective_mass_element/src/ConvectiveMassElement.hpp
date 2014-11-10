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

#ifndef ADOLC
  #error "MoFEM need to be compiled with ADOL-C"
#endif 

namespace MoFEM {

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
    int getRule(int order) { return order-1; };
  };
  
  MyVolumeFE feRhs; ///< cauclate right hand side for tetrahedral elements
  MyVolumeFE& getLoopFeRhs() { return feRhs; } ///< get rhs volume element 
  MyVolumeFE feLhs; //< calculate left hand side for tetrahedral elements
  MyVolumeFE& getLoopFeLhs() { return feLhs; } ///< get lhs volume element

  FieldInterface &mField;
  short int tAg;

  ConvectiveMassElement(
    FieldInterface &m_field,short int tag):
    feRhs(m_field),feLhs(m_field),mField(m_field),tAg(tag) {}

  /** \brief data for calulation het conductivity and heat capacity elements
    * \infroup mofem_forces_and_sources 
    */
  struct BlockData {
    double rho0; ///< reference density
    Range tEts; ///< constatins elements in block set
  }; 
  map<int,BlockData> setOfBlocks; ///< maps block set id with appropiate BlockData

  /** \brief common data used by volume elements
    * \infroup mofem_forces_and_sources 
    */
  struct CommonData {
    vector<ublas::vector<double> > spatialPositionsAtGaussPts;
    vector<ublas::matrix<double> > spatialGardientAtGaussPts;
    vector<ublas::vector<double> > velocityAtGaussPts;
    vector<ublas::matrix<double> > velocityGradientAtGaussPts;
    vector<ublas::vector<double> > meshPositionAtGaussPts;
    vector<ublas::matrix<double> > meshPositionGradientAtGaussPts;
  };
  CommonData commonData;

  /// \brief operator to calulete temeperature gradient at Gauss points
  struct OpGetSpatialAtGaussPts: public TetElementForcesAndSourcesCore::UserDataOperator {

    CommonData &commonData;
    OpGetSpatialAtGaussPts(const string field_name,CommonData &common_data):
      TetElementForcesAndSourcesCore::UserDataOperator(field_name),
      commonData(common_data) {}

    /** \brief operator calulating deformation gradient
      *
      * temerature gradient is calculated multiplying direvatives of shape functions by degrees of freedom
      */
    PetscErrorCode doWork(
      int side,EntityType type,DataForcesAndSurcesCore::EntData &data) {
      PetscFunctionBegin;
      try {

        if(data.getIndices().size()==0) PetscFunctionReturn(0);
 	int nb_dofs = data.getFieldData().size();
	int nb_gauss_pts = data.getN().size1();

	ublas::vector<double>& spatialPositions = data.getFieldData();

	//initialize
	commonData.spatialPositionsAtGaussPts.resize(nb_gauss_pts);
	commonData.spatialGardientAtGaussPts.resize(nb_gauss_pts);
	if(type == MBVERTEX) {
	  for(int gg = 0;gg<nb_gauss_pts;gg++) {
	    commonData.spatialPositionsAtGaussPts[gg].resize(3);
	    commonData.spatialPositionsAtGaussPts[gg].clear();
	    commonData.spatialGardientAtGaussPts[gg].resize(3,3);
	    commonData.spatialGardientAtGaussPts[gg].clear();
	  }
	}	

	for(int gg = 0;gg<nb_gauss_pts;gg++) {
	  ublas::vector<double> N = data.getN(gg,nb_dofs/3);
	  ublas::matrix<double> diffN = trans(data.getDiffN(gg,nb_dofs/3));
	  for(int dd = 0;dd<nb_dofs/3;dd++) {
	    commonData.spatialPositionsAtGaussPts[gg][0] += N[dd]*spatialPositions[3*dd+0];
	    commonData.spatialPositionsAtGaussPts[gg][1] += N[dd]*spatialPositions[3*dd+1];
	    commonData.spatialPositionsAtGaussPts[gg][2] += N[dd]*spatialPositions[3*dd+2];
	    for(int rr = 0;rr<3;rr++) {
	      commonData.spatialGardientAtGaussPts[gg](0,rr) += diffN(rr,dd)*spatialPositions[3*dd+0];
	      commonData.spatialGardientAtGaussPts[gg](1,rr) += diffN(rr,dd)*spatialPositions[3*dd+1];
	      commonData.spatialGardientAtGaussPts[gg](2,rr) += diffN(rr,dd)*spatialPositions[3*dd+2];
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


  /// \brief operator to calulete temeperature gradient at Gauss points
  struct OpGetVelocityGaussPts: public TetElementForcesAndSourcesCore::UserDataOperator {

    CommonData &commonData;
    const EntityType zeroAtType;
    OpGetVelocityGaussPts(const string field_name,CommonData &common_data):
      TetElementForcesAndSourcesCore::UserDataOperator(field_name),
      commonData(common_data),zeroAtType(MBVERTEX) {}

    /** \brief operator calulating deformation gradient
      *
      * temerature gradient is calculated multiplying direvatives of shape functions by degrees of freedom
      */
    PetscErrorCode doWork(
      int side,EntityType type,DataForcesAndSurcesCore::EntData &data) {
      PetscFunctionBegin;
      try {

        if(data.getIndices().size()==0) PetscFunctionReturn(0);
 	int nb_dofs = data.getFieldData().size();
	int nb_gauss_pts = data.getN().size1();

	ublas::vector<double>& spatialVelocities = data.getFieldData();

	//initialize
	commonData.velocityAtGaussPts.resize(nb_gauss_pts);
	commonData.velocityGradientAtGaussPts.resize(nb_gauss_pts);
	if(type == zeroAtType) {
	  for(int gg = 0;gg<nb_gauss_pts;gg++) {
	    commonData.velocityAtGaussPts[gg].resize(3);
	    commonData.velocityAtGaussPts[gg].clear();
	    commonData.velocityGradientAtGaussPts[gg].resize(3,3);
	    commonData.velocityGradientAtGaussPts[gg].clear();
	  }
	}	

	for(int gg = 0;gg<nb_gauss_pts;gg++) {
	  ublas::vector<double> N = data.getN(gg,nb_dofs/3);
	  ublas::matrix<double> diffN = trans(data.getDiffN(gg,nb_dofs/3));
	  for(int dd = 0;dd<nb_dofs/3;dd++) {
	    commonData.velocityAtGaussPts[gg][0] += N[dd]*spatialVelocities[3*dd+0];
	    commonData.velocityAtGaussPts[gg][1] += N[dd]*spatialVelocities[3*dd+1];
	    commonData.velocityAtGaussPts[gg][2] += N[dd]*spatialVelocities[3*dd+2];
	    for(int rr = 0;rr<3;rr++) {
	      commonData.velocityGradientAtGaussPts[gg](0,rr) += diffN(rr,dd)*spatialVelocities[3*dd+0];
	      commonData.velocityGradientAtGaussPts[gg](1,rr) += diffN(rr,dd)*spatialVelocities[3*dd+1];
	      commonData.velocityGradientAtGaussPts[gg](2,rr) += diffN(rr,dd)*spatialVelocities[3*dd+2];
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

  struct OpGetMeshVelocityAtGaussPts: public TetElementForcesAndSourcesCore::UserDataOperator {

    CommonData &commonData;
    OpGetMeshVelocityAtGaussPts(const string field_name,CommonData &common_data):
      TetElementForcesAndSourcesCore::UserDataOperator(field_name),
      commonData(common_data) {}

    /** \brief operator calulating velocity and velocity gradients
      *
      * temerature gradient is calculated multiplying direvatives of shape functions by degrees of freedom
      */
    PetscErrorCode doWork(
      int side,EntityType type,DataForcesAndSurcesCore::EntData &data) {
      PetscFunctionBegin;
      try {

        if(data.getIndices().size()==0) PetscFunctionReturn(0);
 	int nb_dofs = data.getFieldData().size();
	int nb_gauss_pts = data.getN().size1();

	ublas::vector<double>& materialPositions = data.getFieldData();

	//initialize
	commonData.meshPositionAtGaussPts.resize(nb_gauss_pts);
	commonData.meshPositionGradientAtGaussPts.resize(nb_gauss_pts);
	if(type == MBVERTEX) {
	  for(int gg = 0;gg<nb_gauss_pts;gg++) {
	    commonData.meshPositionAtGaussPts[gg].resize(3);
	    commonData.meshPositionGradientAtGaussPts[gg].resize(3,3);
	  }
	}	

	for(int gg = 0;gg<nb_gauss_pts;gg++) {
	  ublas::vector<double> N = data.getN(gg,nb_dofs/3);
	  ublas::matrix<double> diffN = trans(data.getDiffN(gg,nb_dofs/3));
	  for(int dd = 0;dd<nb_dofs/3;dd++) {
	    commonData.meshPositionAtGaussPts[gg][0] += N[dd]*materialPositions[3*dd+0];
	    commonData.meshPositionAtGaussPts[gg][1] += N[dd]*materialPositions[3*dd+1];
	    commonData.meshPositionAtGaussPts[gg][2] += N[dd]*materialPositions[3*dd+2];
	    for(int rr = 0;rr<3;rr++) {
	      commonData.meshPositionGradientAtGaussPts[gg](0,rr) += diffN(rr,dd)*materialPositions[3*dd+0];
	      commonData.meshPositionGradientAtGaussPts[gg](1,rr) += diffN(rr,dd)*materialPositions[3*dd+1];
	      commonData.meshPositionGradientAtGaussPts[gg](2,rr) += diffN(rr,dd)*materialPositions[3*dd+2];
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
      inv_a(0,1) = a(0,3)*a(2,1)-a(0,1)*a(2,2);
      inv_a(0,2) = a(0,1)*a(1,2)-a(0,2)*a(1,1);
      inv_a(1,0) = a(1,2)*a(2,1)-a(1,0)*a(2,2);
      inv_a(1,1) = a(0,0)*a(2,2)-a(0,2)*a(2,0);
      inv_a(1,2) = a(0,2)*a(2,0)-a(0,0)*a(1,2);
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
      ublas::vector<TYPE>& a,ublas::matrix<TYPE>& grad_v,
      ublas::vector<TYPE>& c,ublas::matrix<TYPE>& H,
      ublas::matrix<TYPE>& h,TYPE rho0,
      ublas::vector<TYPE>& dp_dt) {
      PetscFunctionBegin;
  
      PetscErrorCode ierr;
  
      //calulate gradient of deformation
      ublas::matrix<TYPE> invH;
      ierr = iNvert(H,invH); CHKERRQ(ierr);
      ublas::matrix<TYPE> F;
      F = prod(h,invH);
      TYPE detF;
      ierr = dEterminatnt(F,detF); CHKERRQ(ierr);
  
      //calulate current density
      TYPE rho = detF*rho0;
  
      //momentum rate
      dp_dt = rho*(a + prod(grad_v,c));
  
      PetscFunctionReturn(0);
    }
  
    template<typename TYPE> 
    PetscErrorCode calculateFuncUnderIntegral(
      ublas::vector<TYPE>& a,ublas::matrix<TYPE>& grad_v,
      ublas::vector<TYPE>& c,ublas::matrix<TYPE>& H,
      ublas::matrix<TYPE>& h,TYPE rho0,
      ublas::vector<TYPE>& f) {
      PetscFunctionBegin;
  
      PetscErrorCode ierr;
      ublas::vector<TYPE> dp_dt;
      dp_dt.resize(3);
      ierr = calculateMomentumRate(
        a,grad_v,c,H,h,rho0,dp_dt); CHKERRQ(ierr);
      TYPE detH;
      ierr = dEterminatnt(H,detH); CHKERRQ(ierr);
      f = detH*dp_dt;

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
      dAta(data),commonData(common_data) {}
 
    PetscErrorCode doWork(
      int row_side,EntityType row_type,DataForcesAndSurcesCore::EntData &row_data) {
      PetscFunctionBegin;

      PetscErrorCode ierr;
      if(dAta.tEts.find(getMoFEMFEPtr()->get_ent()) == dAta.tEts.end()) {
	PetscFunctionReturn(0);
      }
      if(row_data.getIndices().size()==0) PetscFunctionReturn(0);

      try {

	ublas::vector<double> nf;
	nf.resize(3*row_data.getN().size2(),0);
  
	for(unsigned int gg = 0;gg<row_data.getN().size1();gg++) {

	  ublas::vector<double> a = commonData.velocityAtGaussPts[gg]*getFEMethod()->ts_a;
	  ublas::vector<double> c = commonData.meshPositionAtGaussPts[gg]*getFEMethod()->ts_a;
	  ublas::matrix<double>& g = commonData.velocityGradientAtGaussPts[gg];
	  ublas::matrix<double>& H = commonData.meshPositionGradientAtGaussPts[gg];
	  ublas::matrix<double>& h = commonData.spatialGardientAtGaussPts[gg];

	  double rho0 = dAta.rho0;

	  ublas::vector<double> f;
	  f.resize(3);
	  ierr = calculateFuncUnderIntegral(a,g,c,H,h,rho0,f); CHKERRQ(ierr);
	  double val = getVolume()*getGaussPts()(3,gg);
	  f *= val;

	  for(unsigned int dd = 0;dd<row_data.getN().size2();dd++) {
	    for(int rr = 0;rr<3;rr++) {
	      nf[3*dd+rr] += row_data.getN()(gg,dd)*f[rr];
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

  struct OpMassLhs_dM_dX: public TetElementForcesAndSourcesCore::UserDataOperator,CommonFunctions {

    BlockData &dAta;
    CommonData &commonData;
    int tAg;

    OpMassLhs_dM_dX(const string field_name,BlockData &data,CommonData &common_data,int tag):
      TetElementForcesAndSourcesCore::UserDataOperator(field_name),
      dAta(data),commonData(common_data),tAg(tag) { }

    ublas::vector<adouble> a;
    ublas::matrix<adouble> g;
    ublas::vector<adouble> c;
    ublas::matrix<adouble> H;
    ublas::matrix<adouble> h;

    ublas::vector<double> dX;
    ublas::vector<adouble> a_dX;
    ublas::vector<adouble> a_f;
    double *active_ptr;
    int nb_active_vars;

    virtual PetscErrorCode setPassive(DataForcesAndSurcesCore::EntData &col_data,int gg) {
      PetscFunctionBegin;

      //passive
      H.resize(3,3);
      for(int nn1 = 0;nn1<3;nn1++) {
	for(int nn2 = 0;nn2<3;nn2++) {
	  H(nn1,nn2) = commonData.meshPositionGradientAtGaussPts[gg](nn1,nn2);
	}
      }
      a.resize(3);
      for(int nn = 0;nn<3;nn++) {
	a[nn] = commonData.velocityAtGaussPts[gg][nn]*getFEMethod()->ts_a;
      }
      g.resize(3,3);
      for(int nn1 = 0;nn1<3;nn1++) {
	for(int nn2 = 0;nn2<3;nn2++) {
	  g(nn1,nn2) = commonData.velocityGradientAtGaussPts[gg](nn1,nn2);
	}
      }
      c.resize(3);
      for(int nn = 0;nn<3;nn++) {
	c[nn] = commonData.meshPositionAtGaussPts[gg][nn]*getFEMethod()->ts_a;
      }
      h.resize(3,3);
      for(int nn1 = 0;nn1<3;nn1++) {
	for(int nn2 = 0;nn2<3;nn2++) {
	  h(nn1,nn2) = commonData.spatialGardientAtGaussPts[gg](nn1,nn2);
	}
      }

      PetscFunctionReturn(0);
    }


    virtual PetscErrorCode setActive(DataForcesAndSurcesCore::EntData &col_data,int gg) {
      PetscFunctionBegin;
      
      //active
      dX.resize(col_data.getIndices().size(),0);
      a_dX.resize(dX.size());
      for(unsigned int nn = 0;nn<dX.size();nn++) {
	a_dX[nn] <<= dX[nn];
      }
      for(unsigned int nn1 = 0;nn1<3;nn1++) {
	for(unsigned int nn2 = 0;nn2<3;nn2++) {
	  ublas::matrix<double> diffN = trans(col_data.getDiffN(gg,a_dX.size()));
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
      if(row_data.getIndices().size()==0) PetscFunctionReturn(0);
      if(col_data.getIndices().size()==0) PetscFunctionReturn(0);

      try {

	vector<double*> jac_row_ptr;
	ublas::vector<double> f;
	ublas::matrix<double> jac;
	ublas::matrix<double> k;
	for(unsigned int gg = 0;gg<row_data.getN().size1();gg++) {

	  trace_on(tAg);

	  //set active and passive variables
	  ierr = setPassive(col_data,gg); CHKERRQ(ierr);
	  ierr = setActive(col_data,gg); CHKERRQ(ierr);

	  adouble rho0 = dAta.rho0;
	  a_f.resize(3);
	  ierr = calculateFuncUnderIntegral(a,g,c,H,h,rho0,a_f); CHKERRQ(ierr);

	  //dependant
	  f.resize(3);
	  for(int rr = 0;rr<3;rr++) {
	    a_f[rr] >>= f[rr];
	  }

	  trace_off();

	  if(gg == 0) {
	    jac_row_ptr.resize(a_dX.size());
	    jac.resize(nb_active_vars,a_f.size());
	    for(int nn1 = 0;nn1<nb_active_vars;nn1++) {
	      jac_row_ptr[nn1] = &jac(nn1,0);     
	    }
	    k.resize(row_data.getIndices().size(),col_data.getIndices().size(),0);
	    k.clear();
	  }

	  int r;
	  r = jacobian(
	    tAg,nb_active_vars,a_f.size(),
	    active_ptr,&jac_row_ptr[0]);

	  double val = getVolume()*getGaussPts()(3,gg);
	  for(unsigned int dd1 = 0;dd1<row_data.getN().size2();dd1++) {
	    for(int rr1 = 0;rr1<3;rr1++) {
	      for(unsigned int dd2 = 0;dd2<col_data.getN().size2();dd2++) {
		for(int rr2 = 0;rr2<3;rr2++) {
		  k(3*dd1+rr1,3*dd2+rr2) += val*row_data.getN()(gg,dd1)*jac(3*dd2+rr2,rr1);
		}
	      }
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

  struct OpMassLhs_dM_dx: public OpMassLhs_dM_dX  {

    OpMassLhs_dM_dx(const string field_name,BlockData &data,CommonData &common_data,int tag):
      OpMassLhs_dM_dX(field_name,data,common_data,tag) {}

    ublas::vector<double> dx;
    ublas::vector<adouble> a_dx;

    virtual PetscErrorCode setActive(DataForcesAndSurcesCore::EntData &col_data,int gg) {
      PetscFunctionBegin;
      
      //active
      dx.resize(col_data.getIndices().size(),0);
      a_dx.resize(dx.size());
      for(unsigned int nn = 0;nn<a_dX.size();nn++) {
	a_dx[nn] <<= dx[nn];
      }
      for(unsigned int nn1 = 0;nn1<3;nn1++) {
	for(unsigned int nn2 = 0;nn2<3;nn2++) {
	  ublas::matrix<double> diffN = trans(col_data.getDiffN(gg,a_dx.size()));
	  for(unsigned int dd = 0;dd<a_dx.size();dd++) {
	    h(nn1,nn2) += diffN(nn2,dd)*a_dx[3*dd+nn1];
	  }
	}
      }
      active_ptr = &*dx.data().begin();
      nb_active_vars = dx.size();

      PetscFunctionReturn(0);
    }

  };


  struct OpMassLhs_dM_dc: public OpMassLhs_dM_dX  {

    OpMassLhs_dM_dc(const string field_name,BlockData &data,CommonData &common_data,int tag):
      OpMassLhs_dM_dX(field_name,data,common_data,tag) {}

    ublas::vector<double> dc;
    ublas::vector<adouble> a_dc;

    virtual PetscErrorCode setActive(DataForcesAndSurcesCore::EntData &col_data,int gg) {
      PetscFunctionBegin;

      //active
      dc.resize(col_data.getIndices().size(),0);
      a_dc.resize(dc.size());
      for(unsigned int nn = 0;nn<a_dX.size();nn++) {
	a_dc[nn] <<= dc[nn];
      }
      for(int nn = 0;nn<3;nn++) {
	ublas::vector<double> N = col_data.getN(gg,col_data.getIndices().size()/3);
	for(unsigned int dd = 0;dd<N.size();dd++) {
	  c[nn] += N[dd]*a_dc[3*dd+nn]*getFEMethod()->ts_a;
	}
      }
      active_ptr = &*dc.data().begin();
      nb_active_vars = dc.size();


      PetscFunctionReturn(0);
    }

  };

  struct OpMassLhs_dM_da: public OpMassLhs_dM_dX  {

    OpMassLhs_dM_da(const string field_name,BlockData &data,CommonData &common_data,int tag):
      OpMassLhs_dM_dX(field_name,data,common_data,tag) {}

    ublas::vector<double> dv;
    ublas::vector<adouble> a_dv;

    virtual PetscErrorCode setActive(DataForcesAndSurcesCore::EntData &col_data,int gg) {
      PetscFunctionBegin;

      //active
      dv.resize(col_data.getIndices().size(),0);
      a_dv.resize(dv.size());
      for(unsigned int nn = 0;nn<a_dX.size();nn++) {
	a_dv[nn] <<= dv[nn];
      }
      for(int nn = 0;nn<3;nn++) {
	ublas::vector<double> N = col_data.getN(gg,col_data.getIndices().size()/3);
	for(unsigned int dd = 0;dd<N.size();dd++) {
	  a[nn] += N[dd]*a_dv[3*dd+nn]*getFEMethod()->ts_a;
	}
      }
      active_ptr = &*dv.data().begin();
      nb_active_vars = dv.size();

      PetscFunctionReturn(0);
    }

  };

  struct OpMassLhs_dM_dv: public OpMassLhs_dM_dX  {

    OpMassLhs_dM_dv(const string field_name,BlockData &data,CommonData &common_data,int tag):
      OpMassLhs_dM_dX(field_name,data,common_data,tag) {}

    ublas::vector<double> dv;
    ublas::vector<adouble> a_dv;

    virtual PetscErrorCode setActive(DataForcesAndSurcesCore::EntData &col_data,int gg) {
      PetscFunctionBegin;

      //active
      dv.resize(col_data.getIndices().size(),0);
      a_dv.resize(dv.size());
      for(unsigned int nn = 0;nn<a_dX.size();nn++) {
	a_dv[nn] <<= dv[nn];
      }
      for(unsigned int nn1 = 0;nn1<3;nn1++) {
	for(unsigned int nn2 = 0;nn2<3;nn2++) {
	  ublas::matrix<double> diffN = trans(col_data.getDiffN(gg,a_dv.size()));
	  for(unsigned int dd = 0;dd<a_dv.size();dd++) {
	    g(nn1,nn2) += diffN(nn2,dd)*a_dv[3*dd+nn1];
	  }
	}
      }
      active_ptr = &*dv.data().begin();
      nb_active_vars = dv.size();

      PetscFunctionReturn(0);
    }

  };

  struct OpVelocityRhs: public TetElementForcesAndSourcesCore::UserDataOperator,CommonFunctions {

    BlockData &dAta;
    CommonData &commonData;

    OpVelocityRhs(const string field_name,BlockData &data,CommonData &common_data):
      TetElementForcesAndSourcesCore::UserDataOperator(field_name),
      dAta(data),commonData(common_data) {}
 
    PetscErrorCode doWork(
      int row_side,EntityType row_type,DataForcesAndSurcesCore::EntData &row_data) {
      PetscFunctionBegin;

      PetscErrorCode ierr;
      if(dAta.tEts.find(getMoFEMFEPtr()->get_ent()) == dAta.tEts.end()) {
	PetscFunctionReturn(0);
      }
      if(row_data.getIndices().size()==0) PetscFunctionReturn(0);

      try {

	ublas::vector<double> nf;
	nf.resize(3*row_data.getN().size2(),0);
  
	for(unsigned int gg = 0;gg<row_data.getN().size1();gg++) {

	  ublas::vector<double> v = commonData.velocityAtGaussPts[gg];
	  ublas::vector<double> dot_w = commonData.spatialPositionsAtGaussPts[gg]*getFEMethod()->ts_a;
	  ublas::vector<double> dot_W = commonData.meshPositionAtGaussPts[gg]*getFEMethod()->ts_a;
	  ublas::matrix<double>& H = commonData.meshPositionGradientAtGaussPts[gg];
	  ublas::matrix<double>& h = commonData.spatialGardientAtGaussPts[gg];

	  ublas::vector<double> dot_u;
	  ierr = calulateVelocity(dot_w,dot_W,h,H,dot_u); CHKERRQ(ierr);

	  ublas::vector<double> res = v - dot_u;
	  double val = getVolume()*getGaussPts()(3,gg);
	  double detH;
	  ierr = dEterminatnt(H,detH); CHKERRQ(ierr);
	  val *= detH;

	  for(unsigned int dd = 0;dd<row_data.getN().size2();dd++) {
	    for(int rr = 0;rr<3;rr++) {
	      nf[3*dd+rr] += val*row_data.getN()(gg,dd)*res[rr];
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

  struct OpVelocityLhs_dV_dX: public TetElementForcesAndSourcesCore::UserDataOperator,CommonFunctions {

    BlockData &dAta;
    CommonData &commonData;
    int tAg;

    OpVelocityLhs_dV_dX(const string field_name,BlockData &data,CommonData &common_data,int tag):
      TetElementForcesAndSourcesCore::UserDataOperator(field_name),
      dAta(data),commonData(common_data),tAg(tag) { }

    ublas::vector<adouble> v;
    ublas::vector<adouble> dot_w;
    ublas::vector<adouble> dot_W;
    ublas::matrix<adouble> H;
    ublas::matrix<adouble> h;

    ublas::vector<double> dX;
    ublas::vector<adouble> a_dX;
    double *active_ptr;
    int nb_active_vars;

    virtual PetscErrorCode setPassive(DataForcesAndSurcesCore::EntData &col_data,int gg) {
      PetscFunctionBegin;
      v = commonData.velocityAtGaussPts[gg];
      dot_w = commonData.spatialPositionsAtGaussPts[gg]*getFEMethod()->ts_a;
      dot_W = commonData.meshPositionAtGaussPts[gg]*getFEMethod()->ts_a;
      H = commonData.meshPositionGradientAtGaussPts[gg];
      h = commonData.spatialGardientAtGaussPts[gg];
      PetscFunctionReturn(0);
    }

    virtual PetscErrorCode setActive(DataForcesAndSurcesCore::EntData &col_data,int gg) {
      PetscFunctionBegin;

      //active
      dX.resize(col_data.getIndices().size(),0);
      a_dX.resize(dX.size());
      for(unsigned int nn = 0;nn<dX.size();nn++) {
	a_dX[nn] <<= dX[nn];
      }
      for(unsigned int nn1 = 0;nn1<3;nn1++) {
	for(unsigned int nn2 = 0;nn2<3;nn2++) {
	  ublas::matrix<double> diffN = trans(col_data.getDiffN(gg,a_dX.size()));
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
      if(row_data.getIndices().size()==0) PetscFunctionReturn(0);
      if(col_data.getIndices().size()==0) PetscFunctionReturn(0);

      try {

	vector<double*> jac_row_ptr;
	ublas::matrix<double> jac;
	ublas::matrix<double> k;
	ublas::vector<double> res;
	ublas::vector<adouble> a_res;

	for(unsigned int gg = 0;gg<row_data.getN().size1();gg++) {

	  trace_on(tAg);

	  //set active and passive variables
	  ierr = setPassive(col_data,gg); CHKERRQ(ierr);
	  ierr = setActive(col_data,gg); CHKERRQ(ierr);

	  ublas::vector<adouble> dot_u;
	  ierr = calulateVelocity(dot_w,dot_W,h,H,dot_u); CHKERRQ(ierr);

	  a_res.resize(3);
	  noalias(a_res) = v - dot_u;
	  adouble val = getVolume()*getGaussPts()(3,gg);
	  adouble detH;
	  ierr = dEterminatnt(H,detH); CHKERRQ(ierr);
	  val *= detH;

	  //dependant
	  res.resize(3);
	  for(int rr = 0;rr<3;rr++) {
	    a_res[rr] >>= res[rr];
	  }

	  trace_off();

	  if(gg == 0) {
	    jac_row_ptr.resize(a_dX.size());
	    jac.resize(nb_active_vars,3);
	    for(int nn1 = 0;nn1<nb_active_vars;nn1++) {
	      jac_row_ptr[nn1] = &jac(nn1,0);     
	    }
	    k.resize(row_data.getIndices().size(),col_data.getIndices().size(),0);
	    k.clear();
	  }

	  int r;
	  r = jacobian(
	    tAg,nb_active_vars,3,
	    active_ptr,&jac_row_ptr[0]);

	  { //integrate element stiffnes matrix
	    double val = getVolume()*getGaussPts()(3,gg);
	    for(unsigned int dd1 = 0;dd1<row_data.getN().size2();dd1++) {
	      for(int rr1 = 0;rr1<3;rr1++) {
		for(unsigned int dd2 = 0;dd2<col_data.getN().size2();dd2++) {
		  for(int rr2 = 0;rr2<3;rr2++) {
		    k(3*dd1+rr1,3*dd2+rr2) += val*row_data.getN()(gg,dd1)*jac(3*dd2+rr2,rr1);
		  }
		}
	      }
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

  struct OpVelocityLhs_dV_dx: public OpVelocityLhs_dV_dX {

    OpVelocityLhs_dV_dx(const string field_name,BlockData &data,CommonData &common_data,int tag):
      OpVelocityLhs_dV_dX(field_name,data,common_data,tag) {}

    ublas::vector<adouble> a_dx;
    ublas::vector<double> dx;

    PetscErrorCode setActive(DataForcesAndSurcesCore::EntData &col_data,int gg) {
      PetscFunctionBegin;

      //active
      dx.resize(col_data.getIndices().size(),0);
      a_dx.resize(dX.size());
      for(unsigned int nn = 0;nn<dX.size();nn++) {
	a_dx[nn] <<= dx[nn];
      }
      for(unsigned int nn1 = 0;nn1<3;nn1++) {
	for(unsigned int nn2 = 0;nn2<3;nn2++) {
	  ublas::matrix<double> diffN = trans(col_data.getDiffN(gg,a_dX.size()));
	  for(unsigned int dd = 0;dd<a_dX.size();dd++) {
	    h(nn1,nn2) += diffN(nn2,dd)*a_dx[3*dd+nn1];
	  }
	}
      }
      active_ptr = &*dx.data().begin();
      nb_active_vars = dx.size();

      PetscFunctionReturn(0);
    }

  };

  struct OpVelocityLhs_dV_dot_W: public OpVelocityLhs_dV_dX {

    OpVelocityLhs_dV_dot_W(const string field_name,BlockData &data,CommonData &common_data,int tag):
      OpVelocityLhs_dV_dX(field_name,data,common_data,tag) {}

    PetscErrorCode setActive(DataForcesAndSurcesCore::EntData &col_data,int gg) {
      PetscFunctionBegin;

      //active
      dX.resize(col_data.getIndices().size(),0);
      a_dX.resize(dX.size());
      for(unsigned int nn = 0;nn<dX.size();nn++) {
	a_dX[nn] <<= dX[nn];
      }

      ublas::vector<double> N = col_data.getN(gg,dX.size()/3);
      for(unsigned int nn1 = 0;nn1<3;nn1++) {
	for(unsigned int dd = 0;dd<a_dX.size();dd++) {
	  dot_W[nn1] += N[dd]*a_dX[3*dd+nn1]*getFEMethod()->ts_a;
	}
      }

      active_ptr = &*dX.data().begin();
      nb_active_vars = dX.size();

      PetscFunctionReturn(0);
    }

  };

  struct OpVelocityLhs_dV_dot_w: public OpVelocityLhs_dV_dX {

    OpVelocityLhs_dV_dot_w(const string field_name,BlockData &data,CommonData &common_data,int tag):
      OpVelocityLhs_dV_dX(field_name,data,common_data,tag) {}

    ublas::vector<adouble> a_dx;
    ublas::vector<double> dx;

    PetscErrorCode setActive(DataForcesAndSurcesCore::EntData &col_data,int gg) {
      PetscFunctionBegin;

      //active
      dx.resize(col_data.getIndices().size(),0);
      a_dx.resize(dx.size());
      for(unsigned int nn = 0;nn<dx.size();nn++) {
	a_dx[nn] <<= dx[nn];
      }

      ublas::vector<double> N = col_data.getN(gg,dx.size()/3);
      for(unsigned int nn1 = 0;nn1<3;nn1++) {
	for(unsigned int dd = 0;dd<a_dX.size();dd++) {
	  dot_W[nn1] += N[dd]*a_dX[3*dd+nn1]*getFEMethod()->ts_a;
	}
      }

      active_ptr = &*dx.data().begin();
      nb_active_vars = dx.size();

      PetscFunctionReturn(0);
    }

  };

  PetscErrorCode setBlocks() {
    PetscFunctionBegin;
    ErrorCode rval;
    PetscErrorCode ierr;
  
    for(_IT_CUBITMESHSETS_BY_BCDATA_TYPE_FOR_LOOP_(mField,BLOCKSET|MAT_ELASTICSET,it)) {
      int id = it->get_msId();
      EntityHandle meshset = it->get_meshset();
      rval = mField.get_moab().get_entities_by_type(meshset,MBTET,setOfBlocks[id].tEts,true); CHKERR_PETSC(rval);
      Mat_Elastic mydata;
      ierr = it->get_attribute_data_structure(mydata); CHKERRQ(ierr);
      setOfBlocks[id].rho0 = mydata.data.Density;
    }
    PetscFunctionReturn(0);
  }

  PetscErrorCode addConvectiveMassElement(string element_name,
    string velocity_field_name,
    string spatial_position_field_name,
    string material_position_field_name,bool ale = true) {
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
    if(ale) {
      ierr = mField.modify_finite_element_add_field_row(element_name,material_position_field_name); CHKERRQ(ierr);
      ierr = mField.modify_finite_element_add_field_col(element_name,material_position_field_name); CHKERRQ(ierr);
    }
    ierr = mField.modify_finite_element_add_field_data(element_name,material_position_field_name); CHKERRQ(ierr);

    map<int,BlockData>::iterator sit = setOfBlocks.begin();
    for(;sit!=setOfBlocks.end();sit++) {
      ierr = mField.add_ents_to_finite_element_by_TETs(sit->second.tEts,element_name); CHKERRQ(ierr);
    }

    PetscFunctionReturn(0);
  }

  PetscErrorCode addVelocityElement(string element_name,
    string velocity_field_name,
    string spatial_position_field_name,
    string material_position_field_name,bool ale = true) {
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
    if(ale) {
      ierr = mField.modify_finite_element_add_field_row(element_name,material_position_field_name); CHKERRQ(ierr);
      ierr = mField.modify_finite_element_add_field_col(element_name,material_position_field_name); CHKERRQ(ierr);
      ierr = mField.modify_finite_element_add_field_data(element_name,material_position_field_name); CHKERRQ(ierr);
    }

    map<int,BlockData>::iterator sit = setOfBlocks.begin();
    for(;sit!=setOfBlocks.end();sit++) {
      mField.add_ents_to_finite_element_by_TETs(sit->second.tEts,element_name);
    }

    PetscFunctionReturn(0);
  }

  PetscErrorCode setConvectiveMassOperators(
    string velocity_field_name,
    string spatial_position_field_name,
    string material_position_field_name) {
    PetscFunctionBegin;

    //Rhs
    feRhs.get_op_to_do_Rhs().push_back(new OpGetVelocityGaussPts(velocity_field_name,commonData));
    feRhs.get_op_to_do_Rhs().push_back(new OpGetSpatialAtGaussPts(spatial_position_field_name,commonData));
    feRhs.get_op_to_do_Rhs().push_back(new OpGetMeshVelocityAtGaussPts(material_position_field_name,commonData));
    map<int,BlockData>::iterator sit = setOfBlocks.begin();
    for(;sit!=setOfBlocks.end();sit++) {
      feRhs.get_op_to_do_Rhs().push_back(new OpMassRhs(spatial_position_field_name,sit->second,commonData));
    }

    //Lhs
    feLhs.get_op_to_do_Lhs().push_back(new OpGetVelocityGaussPts(velocity_field_name,commonData));
    feLhs.get_op_to_do_Lhs().push_back(new OpGetSpatialAtGaussPts(spatial_position_field_name,commonData));
    feLhs.get_op_to_do_Lhs().push_back(new OpGetMeshVelocityAtGaussPts(material_position_field_name,commonData));
    sit = setOfBlocks.begin();
    for(;sit!=setOfBlocks.end();sit++) {
    }

    PetscFunctionReturn(0);
  }

  PetscErrorCode setVelocityOperators(
    string velocity_field_name,
    string spatial_position_field_name,
    string material_position_field_name) {
    PetscFunctionBegin;

    //Rhs
    feRhs.get_op_to_do_Rhs().push_back(new OpGetVelocityGaussPts(velocity_field_name,commonData));
    feRhs.get_op_to_do_Rhs().push_back(new OpGetSpatialAtGaussPts(spatial_position_field_name,commonData));
    feRhs.get_op_to_do_Rhs().push_back(new OpGetMeshVelocityAtGaussPts(material_position_field_name,commonData));
    map<int,BlockData>::iterator sit = setOfBlocks.begin();
    for(;sit!=setOfBlocks.end();sit++) {
      feRhs.get_op_to_do_Rhs().push_back(new OpMassRhs(spatial_position_field_name,sit->second,commonData));
    }

    //Lhs
    feLhs.get_op_to_do_Lhs().push_back(new OpGetVelocityGaussPts(velocity_field_name,commonData));
    feLhs.get_op_to_do_Lhs().push_back(new OpGetSpatialAtGaussPts(spatial_position_field_name,commonData));
    feLhs.get_op_to_do_Lhs().push_back(new OpGetMeshVelocityAtGaussPts(material_position_field_name,commonData));
    sit = setOfBlocks.begin();
    for(;sit!=setOfBlocks.end();sit++) {
    }

    PetscFunctionReturn(0);
  }

};

}

#endif //__CONVECTIVE_MASS_ELEMENT_HPP

/***************************************************************************//**
 * \defgroup convective_mass_elem (Convective) Mass Eleement
 * \ingroup mofem_forces_and_sources 
 ******************************************************************************/



