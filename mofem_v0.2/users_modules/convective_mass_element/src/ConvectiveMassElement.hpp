/** \file ConvectiveMassElement.hpp 
 * \brief Operators and data structures for mass and convective mass element
 * \ingroup convective_mass_elem
 *
 */

/* Implementation of convective mass element
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

/** \brief structure grouping operators and data used for calculation of mass (convective) element
  * \ingroup convective_mass_elem
  * \ingroup nonlinear_elastic_elem
  *
  * In order to assemble matrices and right hand vectors, the loops over
  * elements, entities over that elements and finally loop over integration
  * points are executed.
  *
  * Following implementation separate those three celeries of loops and to each
  * loop attach operator.
  *
  */
struct ConvectiveMassElement {

  /// \brief  definition of volume element
  struct MyVolumeFE: public TetElementForcesAndSourcesCore {

    Mat A;
    Vec F;
    bool initV; ///< check if ghost vector used to accumalte Kinetin energy is created

    MyVolumeFE(FieldInterface &_mField): 
      TetElementForcesAndSourcesCore(_mField),A(PETSC_NULL),F(PETSC_NULL),initV(false) {
      meshPositionsFieldName = "NoNE";
    }
    
    /** \brief it is used to calculate nb. of Gauss integration points
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

    Vec V;
    double eNergy;

    PetscErrorCode preProcess() {
      PetscFunctionBegin;
      PetscErrorCode ierr;

      ierr = TetElementForcesAndSourcesCore::preProcess(); CHKERRQ(ierr);

      if(A != PETSC_NULL) {
	ts_B = A;
      }

      if(F != PETSC_NULL) {
	ts_F = F;
      }

      int ghosts[] = { 0 };
      int rank;
      MPI_Comm_rank(mField.get_comm(),&rank);

      switch (ts_ctx) {
	case CTX_TSNONE:
	  if(!initV) {
	    if(rank == 0) {
	      ierr = VecCreateGhost(mField.get_comm(),1,1,1,ghosts,&V); CHKERRQ(ierr);
	    } else {
	      ierr = VecCreateGhost(mField.get_comm(),0,1,1,ghosts,&V); CHKERRQ(ierr);
	    }
	    initV = true;
	  }
	  ierr = VecZeroEntries(V); CHKERRQ(ierr);
	  ierr = VecGhostUpdateBegin(V,INSERT_VALUES,SCATTER_FORWARD); CHKERRQ(ierr);
	  ierr = VecGhostUpdateEnd(V,INSERT_VALUES,SCATTER_FORWARD); CHKERRQ(ierr);
	  break;
	default:
	  break;
      }

      PetscFunctionReturn(0);
    }


    PetscErrorCode postProcess() {
      PetscFunctionBegin;
      PetscErrorCode ierr;

      ierr = TetElementForcesAndSourcesCore::postProcess(); CHKERRQ(ierr);

      double *array;
      switch (ts_ctx) {
	case CTX_TSNONE:
	  ierr = VecAssemblyBegin(V); CHKERRQ(ierr);
	  ierr = VecAssemblyEnd(V); CHKERRQ(ierr);
	  ierr = VecGhostUpdateBegin(V,ADD_VALUES,SCATTER_REVERSE); CHKERRQ(ierr);
	  ierr = VecGhostUpdateEnd(V,ADD_VALUES,SCATTER_REVERSE); CHKERRQ(ierr);
	  ierr = VecGhostUpdateBegin(V,INSERT_VALUES,SCATTER_FORWARD); CHKERRQ(ierr);
	  ierr = VecGhostUpdateEnd(V,INSERT_VALUES,SCATTER_FORWARD); CHKERRQ(ierr);
	  ierr = VecGetArray(V,&array); CHKERRQ(ierr);
	  eNergy = array[0];
	  ierr = VecRestoreArray(V,&array); CHKERRQ(ierr);
	  if(initV) {
	    ierr = VecDestroy(&V); CHKERRQ(ierr);
	    initV = false;
	  }
	  break;
	default:
	  break;
      }

      PetscFunctionReturn(0);
    }


  };

  MyVolumeFE feMassRhs; ///< calculate right hand side for tetrahedral elements
  MyVolumeFE& getLoopFeMassRhs() { return feMassRhs; } ///< get rhs volume element 
  MyVolumeFE feMassLhs; ///< calculate left hand side for tetrahedral elements,i.e. mass element
  MyVolumeFE& getLoopFeMassLhs() { return feMassLhs; } ///< get lhs volume element
  MyVolumeFE feMassAuxLhs; ///< calculate left hand side for tetrahedral elements for Kuu shell matrix
  MyVolumeFE& getLoopFeMassAuxLhs() { return feMassAuxLhs; } ///< get lhs volume element for Kuu shell matrix

  MyVolumeFE feVelRhs; ///< calculate right hand side for tetrahedral elements
  MyVolumeFE& getLoopFeVelRhs() { return feVelRhs; } ///< get rhs volume element 
  MyVolumeFE feVelLhs; ///< calculate left hand side for tetrahedral elements
  MyVolumeFE& getLoopFeVelLhs() { return feVelLhs; } ///< get lhs volume element

  MyVolumeFE feTRhs; ///< calculate right hand side for tetrahedral elements
  MyVolumeFE& getLoopFeTRhs() { return feTRhs; } ///< get rhs volume element 
  MyVolumeFE feTLhs; ///< calculate left hand side for tetrahedral elements
  MyVolumeFE& getLoopFeTLhs() { return feTLhs; } ///< get lhs volume element

  MyVolumeFE feEnergy; ///< calculate kinetic energy
  MyVolumeFE& getLoopFeEnergy() { return feEnergy; } ///< get kinetic energy element

  FieldInterface &mField;
  short int tAg;

  ConvectiveMassElement(
    FieldInterface &m_field,short int tag):
    feMassRhs(m_field),feMassLhs(m_field),feMassAuxLhs(m_field),
    feVelRhs(m_field),feVelLhs(m_field),
    feTRhs(m_field),feTLhs(m_field),
    feEnergy(m_field),
    mField(m_field),tAg(tag) {}

  /** \brief data for calculation inertia forces
    * \ingroup mofem_forces_and_sources 
    */
  struct BlockData {
    double rho0; ///< reference density
    ublas::vector<double> a0; ///< constant acceleration
    Range tEts; ///< elements in block set
  }; 
  map<int,BlockData> setOfBlocks; ///< maps block set id with appropriate BlockData

  /** \brief common data used by volume elements
    * \ingroup mofem_forces_and_sources 
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
    vector<ublas::vector<double> > valT;
    vector<vector<double*> > jacTRowPtr;
    vector<ublas::matrix<double> > jacT;

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

    /** \brief operator calculating deformation gradient
      *
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
      inv_a /= det;
      PetscFunctionReturn(0);
    }
  
  };

  struct OpMassJacobian: public TetElementForcesAndSourcesCore::UserDataOperator,CommonFunctions {

    BlockData &dAta;
    CommonData &commonData;
    int tAg;
    bool jAcobian;
    bool lInear;
    bool fieldDisp;

    OpMassJacobian(const string field_name,BlockData &data,CommonData &common_data,int tag,bool jacobian = true,bool linear = false):
      TetElementForcesAndSourcesCore::UserDataOperator(field_name),
      dAta(data),commonData(common_data),tAg(tag),jAcobian(jacobian),lInear(linear),fieldDisp(false) { }

    ublas::vector<adouble> a,dot_W,dp_dt,a_res;
    ublas::matrix<adouble> h,H,invH,F,g,G;

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
	G.resize(3,3);
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
		if(fieldDisp) {
		  if(nn1==nn2) {
		    h(nn1,nn2) += 1;
		  }
		}
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
	    adouble detH;
	    ierr = dEterminatnt(H,detH); CHKERRQ(ierr);
	    ierr = iNvert(detH,H,invH); CHKERRQ(ierr);
	    noalias(G) = prod(g,invH);
	    double rho0 = dAta.rho0;
	    ublas::vector<double>& a0 = dAta.a0;
	    if(!lInear) {
	      noalias(F) = prod(h,invH);
	      adouble detF;
	      ierr = dEterminatnt(F,detF); CHKERRQ(ierr);
	      //calulate current density
	      adouble rho = rho0*detF;
	      //momentum rate
	      noalias(dp_dt) = rho*(a0 + a + prod(G,dot_W));
	    } else {
	      noalias(dp_dt) = rho0*(a0 + a + prod(G,dot_W));
	    }
	    noalias(a_res) = dp_dt*detH;
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
	  for(int nn1 = 0;nn1<3;nn1++) { //0
	    active[aa++] = (commonData.dataAtGaussPts["DOT_"+commonData.spatialVelocities][gg])[nn1]; 
	  }
	  for(int nn1 = 0;nn1<3;nn1++) { //3
	    for(int nn2 = 0;nn2<3;nn2++) {
	      if(fieldDisp&&nn1 == nn2) {
		active[aa++] = (commonData.gradAtGaussPts[commonData.spatialPositions][gg])(nn1,nn2)+1;
	      } else {
		active[aa++] = (commonData.gradAtGaussPts[commonData.spatialPositions][gg])(nn1,nn2);
	      }
	    }
	  }
	  if(commonData.dataAtGaussPts["DOT_"+commonData.meshPositions].size()>0) {
	    for(int nn1 = 0;nn1<3;nn1++) { //3+9=12
	      for(int nn2 = 0;nn2<3;nn2++) {
	        active[aa++] = (commonData.gradAtGaussPts[commonData.spatialVelocities][gg])(nn1,nn2);
	      }
	    }  
	    for(int nn1 = 0;nn1<3;nn1++) { //3+9+9=21
	      active[aa++] = (commonData.dataAtGaussPts["DOT_"+commonData.meshPositions][gg])[nn1];
	    }
	    for(int nn1 = 0;nn1<3;nn1++) { //3+9+9+3=24
	      for(int nn2 = 0;nn2<3;nn2++) {
		active[aa++] = (commonData.gradAtGaussPts[commonData.meshPositions][gg])(nn1,nn2);
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
	    if(getHoGaussPtsDetJac().size()>0) {
	      val *= getHoGaussPtsDetJac()[gg]; ///< higher order geometry
	    }
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
	    if(getHoGaussPtsDetJac().size()>0) {
	      val *= getHoGaussPtsDetJac()[gg]; ///< higher order geometry
	    }
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
    Range forcesOnlyOnEntities;

    OpMassLhs_dM_dv(
      const string vel_field,const string field_name,BlockData &data,CommonData &common_data,Range *forcesonlyonentities_ptr = NULL):
      TetElementForcesAndSourcesCore::UserDataOperator(vel_field,field_name),
      dAta(data),commonData(common_data) { 
	symm = false;  
	if(forcesonlyonentities_ptr!=NULL) {
	  forcesOnlyOnEntities = *forcesonlyonentities_ptr;
	}
      }

    ublas::matrix<double> k,jac;

    virtual PetscErrorCode getJac(DataForcesAndSurcesCore::EntData &col_data,int gg) {
      PetscFunctionBegin;
      try {
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
      } catch (const std::exception& ex) {
	ostringstream ss;
	ss << "throw in method: " << ex.what() << endl;
	SETERRQ(PETSC_COMM_SELF,1,ss.str().c_str());
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

	  try {
	    ierr = getJac(col_data,gg); CHKERRQ(ierr);
	  } catch (const std::exception& ex) {
	    ostringstream ss;
	    ss << "throw in method: " << ex.what() << endl;
	    SETERRQ(PETSC_COMM_SELF,1,ss.str().c_str());
	  }

	  { //integrate element stiffnes matrix
	    for(int dd1 = 0;dd1<nb_row/3;dd1++) {
	      for(int rr1 = 0;rr1<3;rr1++) {
		for(int dd2 = 0;dd2<nb_col;dd2++) {
		  k(3*dd1+rr1,dd2) += row_data.getN()(gg,dd1)*jac(rr1,dd2);
		}
	      }
	    }
	  }

	}

	if(!forcesOnlyOnEntities.empty()) {
	  ublas::vector<DofIdx> indices = row_data.getIndices();
	  ublas::vector<const FEDofMoFEMEntity*>& dofs = row_data.getFieldDofs();
	  ublas::vector<const FEDofMoFEMEntity*>::iterator dit = dofs.begin();
	  for(int ii = 0;dit!=dofs.end();dit++,ii++) {
	    if(forcesOnlyOnEntities.find((*dit)->get_ent())==forcesOnlyOnEntities.end()) {
	      indices[ii] = -1;
	    }
	  }
	  ierr = MatSetValues(getFEMethod()->ts_B,
	    nb_row,&indices[0],
	    nb_col,&col_data.getIndices()[0],
	    &k(0,0),ADD_VALUES); CHKERRQ(ierr);
	} else {
	  ierr = MatSetValues(getFEMethod()->ts_B,
	    nb_row,&row_data.getIndices()[0],
	    nb_col,&col_data.getIndices()[0],
	    &k(0,0),ADD_VALUES); CHKERRQ(ierr);
	}

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
      try {
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
      } catch (const std::exception& ex) {
	ostringstream ss;
	ss << "throw in method: " << ex.what() << endl;
	SETERRQ(PETSC_COMM_SELF,1,ss.str().c_str());
      }
      PetscFunctionReturn(0);
    }

  };

  struct OpMassLhs_dM_dX: public OpMassLhs_dM_dv  {

    OpMassLhs_dM_dX(const string field_name,const string col_field,BlockData &data,CommonData &common_data):
      OpMassLhs_dM_dv(field_name,col_field,data,common_data) {}

    PetscErrorCode getJac(DataForcesAndSurcesCore::EntData &col_data,int gg) {
      PetscFunctionBegin;
      try {
      int nb_col = col_data.getIndices().size();
      jac.clear();
      //cerr << commonData.jacVel[gg] << endl;
      //cerr << jac << endl;
      ublas::vector<double> N = col_data.getN(gg,nb_col/3);
      for(int dd = 0;dd<nb_col/3;dd++) {
	for(int nn = 0;nn<3;nn++) {
	  jac(0,3*dd+nn) = commonData.jacMass[gg](0,3+9+9+nn)*N(dd)*getFEMethod()->ts_a; 
	  jac(1,3*dd+nn) = commonData.jacMass[gg](1,3+9+9+nn)*N(dd)*getFEMethod()->ts_a; 
	  jac(2,3*dd+nn) = commonData.jacMass[gg](2,3+9+9+nn)*N(dd)*getFEMethod()->ts_a; 
	}
      }
      ublas::matrix<double> diffN = col_data.getDiffN(gg,nb_col/3);
      for(int dd = 0;dd<nb_col/3;dd++) {
	//h00 //h01 //h02
	jac(0,3*dd+0) += commonData.jacMass[gg](0,3+9+9+3+3*0+0)*diffN(dd,0);
	jac(0,3*dd+0) += commonData.jacMass[gg](0,3+9+9+3+3*0+1)*diffN(dd,1);
	jac(0,3*dd+0) += commonData.jacMass[gg](0,3+9+9+3+3*0+2)*diffN(dd,2);
	//h10 //h11 //h12
	jac(1,3*dd+1) += commonData.jacMass[gg](1,3+9+9+3+3*1+0)*diffN(dd,0);
	jac(1,3*dd+1) += commonData.jacMass[gg](1,3+9+9+3+3*1+1)*diffN(dd,1);
	jac(1,3*dd+1) += commonData.jacMass[gg](1,3+9+9+3+3*1+2)*diffN(dd,2);
	//h20 //h21 //h22
	jac(2,3*dd+2) += commonData.jacMass[gg](2,3+9+9+3+3*2+0)*diffN(dd,0);
	jac(2,3*dd+2) += commonData.jacMass[gg](2,3+9+9+3+3*2+1)*diffN(dd,1);
	jac(2,3*dd+2) += commonData.jacMass[gg](2,3+9+9+3+3*2+2)*diffN(dd,2);
      }
      } catch (const std::exception& ex) {
	ostringstream ss;
	ss << "throw in method: " << ex.what() << endl;
	SETERRQ(PETSC_COMM_SELF,1,ss.str().c_str());
      }
      PetscFunctionReturn(0);
    }

  };

  struct OpEnergy: public TetElementForcesAndSourcesCore::UserDataOperator,CommonFunctions {

    BlockData &dAta;
    CommonData &commonData;
    Vec *Vptr;
    bool lInear;

    OpEnergy(const string field_name,BlockData &data,CommonData &common_data,Vec *v_ptr,bool linear = false):
      TetElementForcesAndSourcesCore::UserDataOperator(field_name),
      dAta(data),commonData(common_data),Vptr(v_ptr),lInear(linear) { }

    ublas::matrix<double> h,H,invH,F;
    ublas::vector<double> v;

    PetscErrorCode doWork(
      int row_side,EntityType row_type,DataForcesAndSurcesCore::EntData &row_data) {
      PetscFunctionBegin;

      PetscErrorCode ierr;
      if(row_type != MBVERTEX) {
	PetscFunctionReturn(0);
      }
      if(dAta.tEts.find(getMoFEMFEPtr()->get_ent()) == dAta.tEts.end()) {
	PetscFunctionReturn(0);
      }
      
      try {

	for(unsigned int gg = 0;gg<row_data.getN().size1();gg++) {
	  double val = getVolume()*getGaussPts()(3,gg);
	  if(getHoGaussPtsDetJac().size()>0) {
	    val *= getHoGaussPtsDetJac()[gg]; ///< higher order geometry
	  }
	  double rho0 = dAta.rho0;
	  double rho;
	  if(lInear) {
	    rho = rho0;
	  } else {
	    h.resize(3,3);
	    noalias(h) = (commonData.gradAtGaussPts[commonData.spatialPositions][gg]);
	    if(commonData.dataAtGaussPts["DOT_"+commonData.meshPositions].size()>0) {
	      H.resize(3,3);
	      noalias(H) = (commonData.gradAtGaussPts[commonData.meshPositions][gg]);
	      double detH;
	      ierr = dEterminatnt(H,detH); CHKERRQ(ierr);
	      invH.resize(3,3);
	      ierr = iNvert(detH,H,invH); CHKERRQ(ierr);
	      F.resize(3,3);
	      noalias(F) = prod(h,invH);
	    } else {
	      F.resize(3,3);
	      noalias(F) = h;
	    }
	    double detF;
	    ierr = dEterminatnt(F,detF); CHKERRQ(ierr);
	    rho = detF*rho0;
	  }
	  v.resize(3);
	  noalias(v) = commonData.dataAtGaussPts[commonData.spatialVelocities][gg];
	  double energy = 0.5*rho*inner_prod(v,v);
	  ierr = VecSetValue(*Vptr,0,val*energy,ADD_VALUES); CHKERRQ(ierr);
	}

      } catch (const std::exception& ex) {
	ostringstream ss;
	ss << "throw in method: " << ex.what() << endl;
	SETERRQ(PETSC_COMM_SELF,1,ss.str().c_str());
      }

      PetscFunctionReturn(0);
    }

  };

  struct OpVelocityJacobian: public TetElementForcesAndSourcesCore::UserDataOperator,CommonFunctions {

    BlockData &dAta;
    CommonData &commonData;
    int tAg;
    bool jAcobian,fieldDisp;

    OpVelocityJacobian(const string field_name,BlockData &data,CommonData &common_data,int tag,bool jacobian = true):
      TetElementForcesAndSourcesCore::UserDataOperator(field_name),
      dAta(data),commonData(common_data),tAg(tag),jAcobian(jacobian),fieldDisp(false) { }

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
	      v[nn1] <<= commonData.dataAtGaussPts[commonData.spatialVelocities][gg][nn1];
	      nb_active_vars++;
	    }
	    for(int nn1 = 0;nn1<3;nn1++) { //3
	      dot_w[nn1] <<= commonData.dataAtGaussPts["DOT_"+commonData.spatialPositions][gg][nn1];
	      nb_active_vars++;
	    }
	    if(commonData.dataAtGaussPts["DOT_"+commonData.meshPositions].size()>0) {
	      for(int nn1 = 0;nn1<3;nn1++) { //3+3 = 6
		for(int nn2 = 0;nn2<3;nn2++) {
		  h(nn1,nn2) <<= commonData.gradAtGaussPts[commonData.spatialPositions][gg](nn1,nn2);
		  if(fieldDisp) {
		    if(nn1==nn2) {
		      h(nn1,nn2) += 1;
		    }
		  }
		  nb_active_vars++;
		}
	      }
	      for(int nn1 = 0;nn1<3;nn1++) { //3+3+9
		dot_W[nn1] <<= commonData.dataAtGaussPts["DOT_"+commonData.meshPositions][gg][nn1];
		nb_active_vars++;
	      }
	    }
	    if(commonData.gradAtGaussPts[commonData.meshPositions].size()>0) {
	      for(int nn1 = 0;nn1<3;nn1++) { //3+3+9+3
		for(int nn2 = 0;nn2<3;nn2++) {
		  H(nn1,nn2) <<= commonData.gradAtGaussPts[commonData.meshPositions][gg](nn1,nn2);
		  nb_active_vars++;
		}
	      }
	    }
	    detH = 1;
	    if(commonData.gradAtGaussPts[commonData.meshPositions].size()>0) {
	      ierr = dEterminatnt(H,detH); CHKERRQ(ierr);
	      ierr = iNvert(detH,H,invH); CHKERRQ(ierr);
	      noalias(F) = prod(h,invH);
	    } else {
	      noalias(F) = h;
	    }
	    noalias(dot_u) = dot_w-prod(F,dot_W);
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
		if(fieldDisp&&nn1 == nn2) {
		  active[aa++] = commonData.gradAtGaussPts[commonData.spatialPositions][gg](nn1,nn2)+1;
		} else {
		  active[aa++] = commonData.gradAtGaussPts[commonData.spatialPositions][gg](nn1,nn2);
		}
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
	    if(getHoGaussPtsDetJac().size()>0) {
	      val *= getHoGaussPtsDetJac()[gg]; ///< higher order geometry
	    }
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
	    if(getHoGaussPtsDetJac().size()>0) {
	      val *= getHoGaussPtsDetJac()[gg]; ///< higher order geometry
	    }
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
    bool fieldDisp;

    OpEshelbyDynamicMaterialMomentumJacobian(const string field_name,BlockData &data,CommonData &common_data,int tag,bool jacobian = true):
      TetElementForcesAndSourcesCore::UserDataOperator(field_name),
      dAta(data),commonData(common_data),tAg(tag),jAcobian(jacobian),fieldDisp(false) {}

    ublas::vector<adouble> a,v,a_T;
    ublas::matrix<adouble> g,H,invH,h,F,G;
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

	a.resize(3);
	v.resize(3);
	g.resize(3,3);	
	G.resize(3,3);
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
	commonData.jacTRowPtr.resize(nb_gauss_pts);
	commonData.jacT.resize(nb_gauss_pts);

	int nb_active_vars = 0;
	for(int gg = 0;gg<nb_gauss_pts;gg++) {

	  if(gg == 0) {

	    trace_on(tAg);

	    for(int nn1 = 0;nn1<3;nn1++) { //0
	      a[nn1] <<= commonData.dataAtGaussPts["DOT_"+commonData.spatialVelocities][gg][nn1]; nb_active_vars++;
	    }

	    for(int nn1 = 0;nn1<3;nn1++) { //3
	      v[nn1] <<= commonData.dataAtGaussPts[commonData.spatialVelocities][gg][nn1]; nb_active_vars++;
	    }
	    for(int nn1 = 0;nn1<3;nn1++) { //3+3
	      for(int nn2 = 0;nn2<3;nn2++) {
		g(nn1,nn2) <<= commonData.gradAtGaussPts[commonData.spatialVelocities][gg](nn1,nn2); nb_active_vars++;
	      }
	    }
	    for(int nn1 = 0;nn1<3;nn1++) { //3+3+9
	      for(int nn2 = 0;nn2<3;nn2++) {
		h(nn1,nn2) <<= commonData.gradAtGaussPts[commonData.spatialPositions][gg](nn1,nn2); nb_active_vars++;
		if(fieldDisp) {
		  if(nn1==nn2) {
		    h(nn1,nn2) += 1;
		  }
		}
	      }
	    }
	    if(commonData.gradAtGaussPts[commonData.meshPositions].size()>0) {
	      for(int nn1 = 0;nn1<3;nn1++) { //3+3+9+9
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
	    noalias(G) = prod(g,invH);
	    double rho0 = dAta.rho0;
	    a_T.resize(3);
	    noalias(a_T) = -rho0*(prod(trans(F),a)+prod(trans(G),v))*detH;
	    commonData.valT[gg].resize(3);
	    for(int nn = 0;nn<3;nn++) {
	      a_T[nn] >>= (commonData.valT[gg])[nn];
	    }
	    trace_off();

	  }

	  active.resize(nb_active_vars);
	  int aa = 0;
	    for(int nn1 = 0;nn1<3;nn1++) { //0
	      active[aa++] = commonData.dataAtGaussPts["DOT_"+commonData.spatialVelocities][gg][nn1];
	    }

	    for(int nn1 = 0;nn1<3;nn1++) { //3
	      active[aa++] = commonData.dataAtGaussPts[commonData.spatialVelocities][gg][nn1];
	    }
	    for(int nn1 = 0;nn1<3;nn1++) { //3+3
	      for(int nn2 = 0;nn2<3;nn2++) {
		active[aa++] = commonData.gradAtGaussPts[commonData.spatialVelocities][gg](nn1,nn2); 
	      }
	    }
	    for(int nn1 = 0;nn1<3;nn1++) { //3+3+9
	      for(int nn2 = 0;nn2<3;nn2++) {
		if(fieldDisp&&nn1 == nn2) {
		  active[aa++] = commonData.gradAtGaussPts[commonData.spatialPositions][gg](nn1,nn2)+1; 
		} else {
		  active[aa++] = commonData.gradAtGaussPts[commonData.spatialPositions][gg](nn1,nn2); 
		}
	      }
	    }
	    if(commonData.gradAtGaussPts[commonData.meshPositions].size()>0) {
	      for(int nn1 = 0;nn1<3;nn1++) { //3+3+9+9
		for(int nn2 = 0;nn2<3;nn2++) {
		  active[aa++] = commonData.gradAtGaussPts[commonData.meshPositions][gg](nn1,nn2); 
		}
	      }
	    }

	  if(!jAcobian) {
	    ublas::vector<double>& res = commonData.valT[gg];
	    if(gg>0) {
	      res.resize(3);
	      int r;
	      r = function(tAg,3,nb_active_vars,&active[0],&res[0]);
	      if(r!=3) { // function is locally analytic
		SETERRQ1(PETSC_COMM_SELF,MOFEM_OPERATION_UNSUCCESSFUL,"ADOL-C function evaluation with error r = %d",r);
	      }
	    } 
	    double val = getVolume()*getGaussPts()(3,gg);
	    if(getHoGaussPtsDetJac().size()>0) {
	      val *= getHoGaussPtsDetJac()[gg]; ///< higher order geometry
	    }
	    res *= val;
	  } else {
	    commonData.jacTRowPtr[gg].resize(3);
	    commonData.jacT[gg].resize(3,nb_active_vars);
	    for(int nn1 = 0;nn1<3;nn1++) {
	      (commonData.jacTRowPtr[gg])[nn1] = &(commonData.jacT[gg](nn1,0));     
	    }
	    int r;
	    r = jacobian(
	      tAg,3,nb_active_vars,
	      &active[0],&(commonData.jacTRowPtr[gg])[0]);
	    if(r!=3) {
	      SETERRQ(PETSC_COMM_SELF,MOFEM_OPERATION_UNSUCCESSFUL,"ADOL-C function evaluation with error");
	    }
	    double val = getVolume()*getGaussPts()(3,gg);
	    if(getHoGaussPtsDetJac().size()>0) {
	      val *= getHoGaussPtsDetJac()[gg]; ///< higher order geometry
	    }
	    commonData.jacT[gg] *= val;
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
    Range forcesOnlyOnEntities;

    OpEshelbyDynamicMaterialMomentumRhs(const string field_name,BlockData &data,CommonData &common_data,
      Range *forcesonlyonentities_ptr):
      TetElementForcesAndSourcesCore::UserDataOperator(field_name),
      dAta(data),commonData(common_data) { 
	if(forcesonlyonentities_ptr!=NULL) {
	  forcesOnlyOnEntities = *forcesonlyonentities_ptr;
	}
    }

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
	  ublas::vector<double>& res = commonData.valT[gg];
	  //cerr << res << endl;
	  for(int dd = 0;dd<nb_dofs/3;dd++) {
	    for(int rr = 0;rr<3;rr++) {
	      nf[3*dd+rr] += row_data.getN()(gg,dd)*res[rr];
	    }
	  }
	}

	if(row_data.getIndices().size() > 3*row_data.getN().size2()) {
	  SETERRQ(PETSC_COMM_SELF,1,"data inconsistency");
	}
	if(!forcesOnlyOnEntities.empty()) {
	  ublas::vector<DofIdx> indices = row_data.getIndices();
	  ublas::vector<const FEDofMoFEMEntity*>& dofs = row_data.getFieldDofs();
	  ublas::vector<const FEDofMoFEMEntity*>::iterator dit = dofs.begin();
	  for(int ii = 0;dit!=dofs.end();dit++,ii++) {
	    if(forcesOnlyOnEntities.find((*dit)->get_ent())==forcesOnlyOnEntities.end()) {
	      //cerr << **dit << endl;
	      indices[ii] = -1;
	    }
	  }
	  //cerr << indices << endl;
	  ierr = VecSetValues(getFEMethod()->ts_F,indices.size(),&indices[0],&nf[0],ADD_VALUES); CHKERRQ(ierr);
	} else {
	  ierr = VecSetValues(getFEMethod()->ts_F,row_data.getIndices().size(),
	    &row_data.getIndices()[0],&nf[0],ADD_VALUES); CHKERRQ(ierr);
	}

      } catch (const std::exception& ex) {
	ostringstream ss;
	ss << "throw in method: " << ex.what() << endl;
	SETERRQ(PETSC_COMM_SELF,1,ss.str().c_str());
      }

      PetscFunctionReturn(0);
    }

  };

  struct OpEshelbyDynamicMaterialMomentumLhs_dv: public OpMassLhs_dM_dv {

    OpEshelbyDynamicMaterialMomentumLhs_dv(
      const string vel_field,const string field_name,BlockData &data,CommonData &common_data,
      Range *forcesonlyonentities_ptr):
      OpMassLhs_dM_dv(vel_field,field_name,data,common_data,forcesonlyonentities_ptr) {};


    virtual PetscErrorCode getJac(DataForcesAndSurcesCore::EntData &col_data,int gg) {
      PetscFunctionBegin;
      int nb_col = col_data.getIndices().size();
      jac.clear();
      //cerr << commonData.jacT[gg] << endl;
      //cerr << jac << endl;
      ublas::vector<double> N = col_data.getN(gg,nb_col/3);
      for(int dd = 0;dd<nb_col/3;dd++) {
	for(int nn = 0;nn<3;nn++) {
	  jac(0,3*dd+nn) = commonData.jacT[gg](0,nn)*N(dd)*getFEMethod()->ts_a; 
	  jac(1,3*dd+nn) = commonData.jacT[gg](1,nn)*N(dd)*getFEMethod()->ts_a; 
	  jac(2,3*dd+nn) = commonData.jacT[gg](2,nn)*N(dd)*getFEMethod()->ts_a; 
	}
      }
      for(int dd = 0;dd<nb_col/3;dd++) {
	for(int nn = 0;nn<3;nn++) {
	  jac(0,3*dd+nn) += commonData.jacT[gg](0,3+nn)*N(dd); 
	  jac(1,3*dd+nn) += commonData.jacT[gg](1,3+nn)*N(dd); 
	  jac(2,3*dd+nn) += commonData.jacT[gg](2,3+nn)*N(dd); 
	}
      }
      ublas::matrix<double> diffN = col_data.getDiffN(gg,nb_col/3);
      for(int dd = 0;dd<nb_col/3;dd++) {
	//h00 //h01 //h02
	jac(0,3*dd+0) += commonData.jacT[gg](0,3+3+3*0+0)*diffN(dd,0);
	jac(0,3*dd+0) += commonData.jacT[gg](0,3+3+3*0+1)*diffN(dd,1);
	jac(0,3*dd+0) += commonData.jacT[gg](0,3+3+3*0+2)*diffN(dd,2);
	//h10 //h11 //h12
	jac(1,3*dd+1) += commonData.jacT[gg](1,3+3+3*1+0)*diffN(dd,0);
	jac(1,3*dd+1) += commonData.jacT[gg](1,3+3+3*1+1)*diffN(dd,1);
	jac(1,3*dd+1) += commonData.jacT[gg](1,3+3+3*1+2)*diffN(dd,2);
	//h20 //h21 //h22
	jac(2,3*dd+2) += commonData.jacT[gg](2,3+3+3*2+0)*diffN(dd,0);
	jac(2,3*dd+2) += commonData.jacT[gg](2,3+3+3*2+1)*diffN(dd,1);
	jac(2,3*dd+2) += commonData.jacT[gg](2,3+3+3*2+2)*diffN(dd,2);
      }
      PetscFunctionReturn(0);
    }
  };

  struct OpEshelbyDynamicMaterialMomentumLhs_dx: public OpEshelbyDynamicMaterialMomentumLhs_dv {

    OpEshelbyDynamicMaterialMomentumLhs_dx(
      const string vel_field,const string field_name,BlockData &data,CommonData &common_data,Range *forcesonlyonentities_ptr):
      OpEshelbyDynamicMaterialMomentumLhs_dv(vel_field,field_name,data,common_data,forcesonlyonentities_ptr) {};

    virtual PetscErrorCode getJac(DataForcesAndSurcesCore::EntData &col_data,int gg) {
      PetscFunctionBegin;
      int nb_col = col_data.getIndices().size();
      jac.clear();
      //cerr << commonData.jacT[gg] << endl;
      //cerr << jac << endl;
      ublas::matrix<double> diffN = col_data.getDiffN(gg,nb_col/3);
      for(int dd = 0;dd<nb_col/3;dd++) {
	//h00 //h01 //h02
	jac(0,3*dd+0) += commonData.jacT[gg](0,3+3+9+3*0+0)*diffN(dd,0);
	jac(0,3*dd+0) += commonData.jacT[gg](0,3+3+9+3*0+1)*diffN(dd,1);
	jac(0,3*dd+0) += commonData.jacT[gg](0,3+3+9+3*0+2)*diffN(dd,2);
	//h10 //h11 //h12
	jac(1,3*dd+1) += commonData.jacT[gg](1,3+3+9+3*1+0)*diffN(dd,0);
	jac(1,3*dd+1) += commonData.jacT[gg](1,3+3+9+3*1+1)*diffN(dd,1);
	jac(1,3*dd+1) += commonData.jacT[gg](1,3+3+9+3*1+2)*diffN(dd,2);
	//h20 //h21 //h22
	jac(2,3*dd+2) += commonData.jacT[gg](2,3+3+9+3*2+0)*diffN(dd,0);
	jac(2,3*dd+2) += commonData.jacT[gg](2,3+3+9+3*2+1)*diffN(dd,1);
	jac(2,3*dd+2) += commonData.jacT[gg](2,3+3+9+3*2+2)*diffN(dd,2);
      }
      PetscFunctionReturn(0);
    }

  };

  struct OpEshelbyDynamicMaterialMomentumLhs_dX: public OpEshelbyDynamicMaterialMomentumLhs_dv {

    OpEshelbyDynamicMaterialMomentumLhs_dX(
      const string vel_field,const string field_name,BlockData &data,CommonData &common_data,Range *forcesonlyonentities_ptr):
      OpEshelbyDynamicMaterialMomentumLhs_dv(vel_field,field_name,data,common_data,forcesonlyonentities_ptr) {};

    virtual PetscErrorCode getJac(DataForcesAndSurcesCore::EntData &col_data,int gg) {
      PetscFunctionBegin;
      int nb_col = col_data.getIndices().size();
      jac.clear();
      //cerr << commonData.jacT[gg] << endl;
      //cerr << jac << endl;
      ublas::matrix<double> diffN = col_data.getDiffN(gg,nb_col/3);
      for(int dd = 0;dd<nb_col/3;dd++) {
	//h00 //h01 //h02
	jac(0,3*dd+0) += commonData.jacT[gg](0,3+3+9+9+3*0+0)*diffN(dd,0);
	jac(0,3*dd+0) += commonData.jacT[gg](0,3+3+9+9+3*0+1)*diffN(dd,1);
	jac(0,3*dd+0) += commonData.jacT[gg](0,3+3+9+9+3*0+2)*diffN(dd,2);
	//h10 //h11 //h12
	jac(1,3*dd+1) += commonData.jacT[gg](1,3+3+9+9+3*1+0)*diffN(dd,0);
	jac(1,3*dd+1) += commonData.jacT[gg](1,3+3+9+9+3*1+1)*diffN(dd,1);
	jac(1,3*dd+1) += commonData.jacT[gg](1,3+3+9+9+3*1+2)*diffN(dd,2);
	//h20 //h21 //h22
	jac(2,3*dd+2) += commonData.jacT[gg](2,3+3+9+9+3*2+0)*diffN(dd,0);
	jac(2,3*dd+2) += commonData.jacT[gg](2,3+3+9+9+3*2+1)*diffN(dd,1);
	jac(2,3*dd+2) += commonData.jacT[gg](2,3+3+9+9+3*2+2)*diffN(dd,2);
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

      //ierr = mField.set_other_local_ghost_vector(problemPtr,velocityField,"DOT_"+velocityField,COL,ts_u_t,INSERT_VALUES,SCATTER_REVERSE); CHKERRQ(ierr);
      //ierr = mField.set_other_local_ghost_vector(problemPtr,spatialPositionField,"DOT_"+spatialPositionField,COL,ts_u_t,INSERT_VALUES,SCATTER_REVERSE); CHKERRQ(ierr);

      //FIXME: This global scattering because Kuu problem and Dynamic problem
      //not share partitions. Both problem should use the same partitioning to
      //resolve this problem.
      ierr = mField.set_global_ghost_vector(problemPtr,COL,ts_u,INSERT_VALUES,SCATTER_REVERSE); CHKERRQ(ierr);
      ierr = mField.set_other_global_ghost_vector(problemPtr,velocityField,"DOT_"+velocityField,COL,ts_u_t,INSERT_VALUES,SCATTER_REVERSE); CHKERRQ(ierr);
      ierr = mField.set_other_global_ghost_vector(problemPtr,spatialPositionField,"DOT_"+spatialPositionField,COL,ts_u_t,INSERT_VALUES,SCATTER_REVERSE); CHKERRQ(ierr);

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
    bool ale = false,
    BitRefLevel bit = BitRefLevel(),Range *intersected = NULL) {
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
    if(intersected!=NULL) {
      if(tets.empty()) {
	tets = *intersected;
      } else {
	tets = intersect(*intersected,tets);
      }
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
    bool ale = false,bool linear = false) {
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
      } else {
	feMassRhs.meshPositionsFieldName = material_position_field_name;
      }
    }
    map<int,BlockData>::iterator sit = setOfBlocks.begin();
    for(;sit!=setOfBlocks.end();sit++) {
      feMassRhs.get_op_to_do_Rhs().push_back(new OpMassJacobian(spatial_position_field_name,sit->second,commonData,tAg,false,linear));
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
      } else {
	feMassLhs.meshPositionsFieldName = material_position_field_name;
      }
    }
    sit = setOfBlocks.begin();
    for(;sit!=setOfBlocks.end();sit++) {
      feMassLhs.get_op_to_do_Rhs().push_back(new OpMassJacobian(spatial_position_field_name,sit->second,commonData,tAg,true,linear));
      feMassLhs.get_op_to_do_Lhs().push_back(new OpMassLhs_dM_dv(spatial_position_field_name,velocity_field_name,sit->second,commonData));
      feMassLhs.get_op_to_do_Lhs().push_back(new OpMassLhs_dM_dx(spatial_position_field_name,spatial_position_field_name,sit->second,commonData));
      if(mField.check_field(material_position_field_name)) {
	if(ale) {
	  feMassLhs.get_op_to_do_Lhs().push_back(new OpMassLhs_dM_dX(spatial_position_field_name,material_position_field_name,sit->second,commonData));
	} else {
	  feMassLhs.meshPositionsFieldName = material_position_field_name;
	}
      }
    }

    //Energy
    feEnergy.get_op_to_do_Rhs().push_back(new OpGetCommonDataAtGaussPts(velocity_field_name,commonData));
    feEnergy.get_op_to_do_Rhs().push_back(new OpGetCommonDataAtGaussPts(spatial_position_field_name,commonData));
    if(mField.check_field(material_position_field_name)) {
      feEnergy.get_op_to_do_Rhs().push_back(new OpGetCommonDataAtGaussPts(material_position_field_name,commonData));
      feEnergy.meshPositionsFieldName = material_position_field_name;
    }
    sit = setOfBlocks.begin();
    for(;sit!=setOfBlocks.end();sit++) {
      feEnergy.get_op_to_do_Rhs().push_back(new OpEnergy(spatial_position_field_name,sit->second,commonData,&feEnergy.V,linear));
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
    if(mField.check_field(material_position_field_name)) {
      feVelRhs.get_op_to_do_Rhs().push_back(new OpGetCommonDataAtGaussPts("DOT_"+spatial_position_field_name,commonData));
      feVelRhs.get_op_to_do_Rhs().push_back(new OpGetCommonDataAtGaussPts(material_position_field_name,commonData));
      if(ale) {
	feVelRhs.get_op_to_do_Rhs().push_back(new OpGetCommonDataAtGaussPts("DOT_"+material_position_field_name,commonData));
      } else {
	feVelRhs.meshPositionsFieldName = material_position_field_name;
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
    if(mField.check_field(material_position_field_name)) {
      feVelLhs.get_op_to_do_Rhs().push_back(new OpGetCommonDataAtGaussPts("DOT_"+spatial_position_field_name,commonData));
      feVelLhs.get_op_to_do_Rhs().push_back(new OpGetCommonDataAtGaussPts(material_position_field_name,commonData));
      if(ale) {
	feVelLhs.get_op_to_do_Rhs().push_back(new OpGetCommonDataAtGaussPts("DOT_"+material_position_field_name,commonData));
      } else {
	feVelLhs.meshPositionsFieldName = material_position_field_name;
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
	} else {
	  feVelLhs.meshPositionsFieldName = material_position_field_name;
	}
      }
    }

    PetscFunctionReturn(0);
  }


  PetscErrorCode setKinematicEshelbyOperators(
    string velocity_field_name,
    string spatial_position_field_name,
    string material_position_field_name = "MESH_NODE_POSITIONS",
    Range *forces_on_entities_ptr = NULL) {
    PetscFunctionBegin;

    commonData.spatialPositions = spatial_position_field_name;
    commonData.meshPositions = material_position_field_name;
    commonData.spatialVelocities = velocity_field_name;

    //Rhs
    feTRhs.get_op_to_do_Rhs().push_back(new OpGetCommonDataAtGaussPts(velocity_field_name,commonData));
    feTRhs.get_op_to_do_Rhs().push_back(new OpGetCommonDataAtGaussPts(spatial_position_field_name,commonData));
    feTRhs.get_op_to_do_Rhs().push_back(new OpGetCommonDataAtGaussPts(material_position_field_name,commonData));
    feTRhs.get_op_to_do_Rhs().push_back(new OpGetCommonDataAtGaussPts("DOT_"+velocity_field_name,commonData));

    map<int,BlockData>::iterator sit = setOfBlocks.begin();
    for(;sit!=setOfBlocks.end();sit++) {
      feTRhs.get_op_to_do_Rhs().push_back(new OpEshelbyDynamicMaterialMomentumJacobian(material_position_field_name,sit->second,commonData,tAg,false));
      feTRhs.get_op_to_do_Rhs().push_back(
	new OpEshelbyDynamicMaterialMomentumRhs(material_position_field_name,sit->second,commonData,forces_on_entities_ptr));
    }

    //Lhs
    feTLhs.get_op_to_do_Rhs().push_back(new OpGetCommonDataAtGaussPts(velocity_field_name,commonData));
    feTLhs.get_op_to_do_Rhs().push_back(new OpGetCommonDataAtGaussPts(spatial_position_field_name,commonData));
    feTLhs.get_op_to_do_Rhs().push_back(new OpGetCommonDataAtGaussPts("DOT_"+velocity_field_name,commonData));
    if(mField.check_field(material_position_field_name)) {
      feTLhs.get_op_to_do_Rhs().push_back(new OpGetCommonDataAtGaussPts(material_position_field_name,commonData));
    }
    sit = setOfBlocks.begin();
    for(;sit!=setOfBlocks.end();sit++) {
      feTLhs.get_op_to_do_Rhs().push_back(new OpEshelbyDynamicMaterialMomentumJacobian(material_position_field_name,sit->second,commonData,tAg));
      feTLhs.get_op_to_do_Lhs().push_back(
	new OpEshelbyDynamicMaterialMomentumLhs_dv(material_position_field_name,velocity_field_name,sit->second,commonData,forces_on_entities_ptr));
      feTLhs.get_op_to_do_Lhs().push_back(
	new OpEshelbyDynamicMaterialMomentumLhs_dx(material_position_field_name,spatial_position_field_name,sit->second,commonData,forces_on_entities_ptr));
      feTLhs.get_op_to_do_Lhs().push_back(
	new OpEshelbyDynamicMaterialMomentumLhs_dX(material_position_field_name,material_position_field_name,sit->second,commonData,forces_on_entities_ptr));
    }

    PetscFunctionReturn(0);
  }

  PetscErrorCode setShellMatrixMassOperators(
    string velocity_field_name,
    string spatial_position_field_name,
    string material_position_field_name = "MESH_NODE_POSITIONS",
    bool linear = false) {
    PetscFunctionBegin;

    commonData.spatialPositions = spatial_position_field_name;
    commonData.meshPositions = material_position_field_name;
    commonData.spatialVelocities = velocity_field_name;

    //Rhs
    feMassRhs.get_op_to_do_Rhs().push_back(new OpGetCommonDataAtGaussPts(velocity_field_name,commonData));
    feMassRhs.get_op_to_do_Rhs().push_back(new OpGetCommonDataAtGaussPts(spatial_position_field_name,commonData));
    feMassRhs.get_op_to_do_Rhs().push_back(new OpGetCommonDataAtGaussPts("DOT_"+velocity_field_name,commonData));
    if(mField.check_field(material_position_field_name)) {
      feMassRhs.get_op_to_do_Rhs().push_back(new OpGetCommonDataAtGaussPts(material_position_field_name,commonData));
      feMassRhs.meshPositionsFieldName = material_position_field_name;
    }
    map<int,BlockData>::iterator sit = setOfBlocks.begin();
    for(;sit!=setOfBlocks.end();sit++) {
      feMassRhs.get_op_to_do_Rhs().push_back(new OpMassJacobian(spatial_position_field_name,sit->second,commonData,tAg,false,linear));
      feMassRhs.get_op_to_do_Rhs().push_back(new OpMassRhs(spatial_position_field_name,sit->second,commonData));
    }

    //Lhs
    feMassLhs.get_op_to_do_Rhs().push_back(new OpGetCommonDataAtGaussPts(velocity_field_name,commonData));
    feMassLhs.get_op_to_do_Rhs().push_back(new OpGetCommonDataAtGaussPts(spatial_position_field_name,commonData));
    feMassLhs.get_op_to_do_Rhs().push_back(new OpGetCommonDataAtGaussPts("DOT_"+velocity_field_name,commonData));
    if(mField.check_field(material_position_field_name)) {
      feMassLhs.get_op_to_do_Rhs().push_back(new OpGetCommonDataAtGaussPts(material_position_field_name,commonData));
      feMassLhs.meshPositionsFieldName = material_position_field_name;
    }
    sit = setOfBlocks.begin();
    for(;sit!=setOfBlocks.end();sit++) {
      feMassLhs.get_op_to_do_Rhs().push_back(new OpMassJacobian(spatial_position_field_name,sit->second,commonData,tAg,true,linear));
      feMassLhs.get_op_to_do_Lhs().push_back(new OpMassLhs_dM_dv(spatial_position_field_name,spatial_position_field_name,sit->second,commonData));
      if(mField.check_field(material_position_field_name)) {
	feMassLhs.meshPositionsFieldName = material_position_field_name;
      }
    }

    //Aux Lhs
    feMassAuxLhs.get_op_to_do_Rhs().push_back(new OpGetCommonDataAtGaussPts(velocity_field_name,commonData));
    feMassAuxLhs.get_op_to_do_Rhs().push_back(new OpGetCommonDataAtGaussPts(spatial_position_field_name,commonData));
    feMassAuxLhs.get_op_to_do_Rhs().push_back(new OpGetCommonDataAtGaussPts("DOT_"+velocity_field_name,commonData));
    if(mField.check_field(material_position_field_name)) {
      feMassAuxLhs.get_op_to_do_Rhs().push_back(new OpGetCommonDataAtGaussPts(material_position_field_name,commonData));
      feMassAuxLhs.meshPositionsFieldName = material_position_field_name;
    }
    sit = setOfBlocks.begin();
    for(;sit!=setOfBlocks.end();sit++) {
      feMassAuxLhs.get_op_to_do_Rhs().push_back(new OpMassJacobian(spatial_position_field_name,sit->second,commonData,tAg,true,linear));
      feMassAuxLhs.get_op_to_do_Lhs().push_back(new OpMassLhs_dM_dx(spatial_position_field_name,spatial_position_field_name,sit->second,commonData));
      if(mField.check_field(material_position_field_name)) {
	feMassAuxLhs.meshPositionsFieldName = material_position_field_name;
      }
    }

    //Energy E=0.5*rho*v*v
    feEnergy.get_op_to_do_Rhs().push_back(new OpGetCommonDataAtGaussPts(velocity_field_name,commonData));
    feEnergy.get_op_to_do_Rhs().push_back(new OpGetCommonDataAtGaussPts(spatial_position_field_name,commonData));
    if(mField.check_field(material_position_field_name)) {
      feEnergy.get_op_to_do_Rhs().push_back(new OpGetCommonDataAtGaussPts(material_position_field_name,commonData));
      feEnergy.meshPositionsFieldName = material_position_field_name;
    }
    sit = setOfBlocks.begin();
    for(;sit!=setOfBlocks.end();sit++) {
      feEnergy.get_op_to_do_Rhs().push_back(new OpEnergy(spatial_position_field_name,sit->second,commonData,&feEnergy.V,linear));
    }

    PetscFunctionReturn(0);
  }


  struct MatShellCtx {

    Mat K,M;
    VecScatter scatterU,scatterV;
    double ts_a;//,scale;

    bool iNitialized;
    MatShellCtx(): iNitialized(false) {}
    virtual ~MatShellCtx() {
      if(iNitialized) {
	PetscErrorCode ierr;
	ierr = dEstroy(); CHKERRABORT(PETSC_COMM_WORLD,ierr);
      }
    }

    Mat barK;
    Vec u,v,Ku,Mv;
    PetscErrorCode iNit() {
      PetscFunctionBegin;
      if(!iNitialized) {
	PetscErrorCode ierr;
	#if PETSC_VERSION_GE(3,5,3) 
	ierr = MatCreateVecs(K,&u,&Ku); CHKERRQ(ierr);
	ierr = MatCreateVecs(M,&v,&Mv); CHKERRQ(ierr);
	#else 
	ierr = MatGetVecs(K,&u,&Ku); CHKERRQ(ierr);
	ierr = MatGetVecs(M,&v,&Mv); CHKERRQ(ierr);
	#endif
	ierr = MatDuplicate(K,MAT_SHARE_NONZERO_PATTERN,&barK); CHKERRQ(ierr);
	iNitialized = true;
      }
      PetscFunctionReturn(0);
    }

    PetscErrorCode dEstroy() {
      PetscFunctionBegin;
      if(iNitialized) {
	PetscErrorCode ierr;
	ierr = VecDestroy(&u); CHKERRQ(ierr);
	ierr = VecDestroy(&Ku); CHKERRQ(ierr);
	ierr = VecDestroy(&v); CHKERRQ(ierr);
	ierr = VecDestroy(&Mv); CHKERRQ(ierr);
	ierr = MatDestroy(&barK); CHKERRQ(ierr);
	iNitialized = false;
      }
      PetscFunctionReturn(0);
    }

    friend PetscErrorCode MultOpA(Mat A,Vec x,Vec f);
    friend PetscErrorCode ZeroEntriesOp(Mat A);
  
  };

  /** \brief Mult operator for shell matrix
    *
    * \f[
    \left[
    \begin{array}{cc}
    \mathbf{M} & \mathbf{K} \\
    \mathbf{I} & -\mathbf{I}a 
    \end{array}
    \right]
    \left[
    \begin{array}{c}
    \mathbf{v} \\
    \mathbf{u}
    \end{array}
    \right] =
    \left[
    \begin{array}{c}
    \mathbf{r}_u \\
    \mathbf{r}_v
    \end{array}
    \right]
    * \f]
    * 
    */
  static PetscErrorCode MultOpA(Mat A,Vec x,Vec f) {
    PetscFunctionBegin;
    PetscErrorCode ierr;
    void *void_ctx;
    ierr = MatShellGetContext(A,&void_ctx); CHKERRQ(ierr);
    MatShellCtx *ctx = (MatShellCtx*)void_ctx;
    if(!ctx->iNitialized) {
      ierr = ctx->iNit(); CHKERRQ(ierr);
    }
    ierr = VecZeroEntries(f); CHKERRQ(ierr);
    //Mult Ku
    ierr = VecScatterBegin(ctx->scatterU,x,ctx->u,INSERT_VALUES,SCATTER_FORWARD); CHKERRQ(ierr);
    ierr = VecScatterEnd(ctx->scatterU,x,ctx->u,INSERT_VALUES,SCATTER_FORWARD); CHKERRQ(ierr);
    ierr = MatMult(ctx->K,ctx->u,ctx->Ku); CHKERRQ(ierr);
    ierr = VecScatterBegin(ctx->scatterU,ctx->Ku,f,INSERT_VALUES,SCATTER_REVERSE); CHKERRQ(ierr);
    ierr = VecScatterEnd(ctx->scatterU,ctx->Ku,f,INSERT_VALUES,SCATTER_REVERSE); CHKERRQ(ierr);
    //Mult Mv
    ierr = VecScatterBegin(ctx->scatterV,x,ctx->v,INSERT_VALUES,SCATTER_FORWARD); CHKERRQ(ierr);
    ierr = VecScatterEnd(ctx->scatterV,x,ctx->v,INSERT_VALUES,SCATTER_FORWARD); CHKERRQ(ierr);
    ierr = MatMult(ctx->M,ctx->v,ctx->Mv); CHKERRQ(ierr);
    ierr = VecScatterBegin(ctx->scatterU,ctx->Mv,f,ADD_VALUES,SCATTER_REVERSE); CHKERRQ(ierr);
    ierr = VecScatterEnd(ctx->scatterU,ctx->Mv,f,ADD_VALUES,SCATTER_REVERSE); CHKERRQ(ierr);
    //Velocities
    ierr = VecAXPY(ctx->v,-ctx->ts_a,ctx->u); CHKERRQ(ierr);
    //ierr = VecScale(ctx->v,ctx->scale); CHKERRQ(ierr);
    ierr = VecScatterBegin(ctx->scatterV,ctx->v,f,INSERT_VALUES,SCATTER_REVERSE); CHKERRQ(ierr);
    ierr = VecScatterEnd(ctx->scatterV,ctx->v,f,INSERT_VALUES,SCATTER_REVERSE); CHKERRQ(ierr);
    //Assemble
    ierr = VecAssemblyBegin(f); CHKERRQ(ierr);
    ierr = VecAssemblyEnd(f); CHKERRQ(ierr);
    PetscFunctionReturn(0);
  }

  static PetscErrorCode ZeroEntriesOp(Mat A) {
    PetscFunctionBegin;
    PetscErrorCode ierr;
    void *void_ctx;
    ierr = MatShellGetContext(A,&void_ctx); CHKERRQ(ierr);
    MatShellCtx *ctx = (MatShellCtx*)void_ctx;
    ierr = MatZeroEntries(ctx->K); CHKERRQ(ierr);
    ierr = MatZeroEntries(ctx->M); CHKERRQ(ierr);
    PetscFunctionReturn(0);
  }

  struct PCShellCtx {

    Mat shellMat;
    bool initPC; ///< check if PC is initialized

    PCShellCtx(Mat shell_mat):
      shellMat(shell_mat),initPC(false) {
    }
 
    PC pC;
  
    PetscErrorCode iNit() {
      PetscFunctionBegin;
      PetscErrorCode ierr;
      if(!initPC) {
	MPI_Comm comm;
	ierr = PetscObjectGetComm((PetscObject)shellMat,&comm); CHKERRQ(ierr);
	ierr = PCCreate(comm,&pC); CHKERRQ(ierr);
	initPC = true;
      }
      PetscFunctionReturn(0);
    }

    PetscErrorCode dEstroy() {
      PetscFunctionBegin;
      PetscErrorCode ierr;
      if(initPC) {
	ierr = PCDestroy(&pC); CHKERRQ(ierr);
	initPC = false;
      }
      PetscFunctionReturn(0);
    }

    friend PetscErrorCode PCShellSetUpOp(PC pc);
    friend PetscErrorCode PCShellDestroy(PC pc);
    friend PetscErrorCode PCShellApplyOp(PC pc,Vec f,Vec x);

  };
  
  static PetscErrorCode PCShellSetUpOp(PC pc) {
    PetscFunctionBegin;
    PetscErrorCode ierr;
    void *void_ctx;
    ierr = PCShellGetContext(pc,&void_ctx); CHKERRQ(ierr);
    PCShellCtx *ctx = (PCShellCtx*)void_ctx;
    ierr = ctx->iNit(); CHKERRQ(ierr);
    MatShellCtx *shell_mat_ctx;
    ierr = MatShellGetContext(ctx->shellMat,&shell_mat_ctx); CHKERRQ(ierr);
    ierr = PCSetFromOptions(ctx->pC); CHKERRQ(ierr);
    ierr = PCSetOperators(ctx->pC,shell_mat_ctx->barK,shell_mat_ctx->barK); CHKERRQ(ierr);
    ierr = PCSetUp(ctx->pC); CHKERRQ(ierr);
    PetscFunctionReturn(0);
  }

  static PetscErrorCode PCShellDestroy(PC pc) {
    PetscFunctionBegin;
    PetscErrorCode ierr;
    void *void_ctx;
    ierr = PCShellGetContext(pc,&void_ctx); CHKERRQ(ierr);
    PCShellCtx *ctx = (PCShellCtx*)void_ctx;
    ierr = ctx->dEstroy(); CHKERRQ(ierr);
    PetscFunctionReturn(0);
  }

  /** \brief apply pre-conditioner for shell matrix 
    *
    * \f[
    \left[
    \begin{array}{cc}
    \mathbf{M} & \mathbf{K} \\
    \mathbf{I} & -\mathbf{I}a 
    \end{array}
    \right]
    \left[
    \begin{array}{c}
    \mathbf{v} \\
    \mathbf{u}
    \end{array}
    \right] =
    \left[
    \begin{array}{c}
    \mathbf{r}_u \\
    \mathbf{r}_v
    \end{array}
    \right]
    * \f]
    *
    * where \f$\mathbf{v} = \mathbf{r}_v + a\mathbf{u}\f$ and \f$\mathbf{u}=(a\mathbf{M}+\mathbf{K})^{-1}(\mathbf{r}_u - \mathbf{M}\mathbf{r}_v\f$.
    *
    */
  static PetscErrorCode PCShellApplyOp(PC pc,Vec f,Vec x) {
    PetscFunctionBegin;
    PetscErrorCode ierr;
    void *void_ctx;
    ierr = PCShellGetContext(pc,&void_ctx); CHKERRQ(ierr);
    PCShellCtx *ctx = (PCShellCtx*)void_ctx;
    MatShellCtx *shell_mat_ctx;
    ierr = MatShellGetContext(ctx->shellMat,&shell_mat_ctx); CHKERRQ(ierr);
    //forward
    ierr = VecScatterBegin(shell_mat_ctx->scatterU,f,shell_mat_ctx->Ku,INSERT_VALUES,SCATTER_FORWARD); CHKERRQ(ierr);
    ierr = VecScatterEnd(shell_mat_ctx->scatterU,f,shell_mat_ctx->Ku,INSERT_VALUES,SCATTER_FORWARD); CHKERRQ(ierr);
    ierr = VecScatterBegin(shell_mat_ctx->scatterV,f,shell_mat_ctx->v,INSERT_VALUES,SCATTER_FORWARD); CHKERRQ(ierr);
    ierr = VecScatterEnd(shell_mat_ctx->scatterV,f,shell_mat_ctx->v,INSERT_VALUES,SCATTER_FORWARD); CHKERRQ(ierr);
    //ierr = VecScale(shell_mat_ctx->v,1/shell_mat_ctx->scale); CHKERRQ(ierr);
    //apply pre-conditioner and calculate u
    ierr = MatMult(shell_mat_ctx->M,shell_mat_ctx->v,shell_mat_ctx->Mv); CHKERRQ(ierr); // Mrv
    ierr = VecAXPY(shell_mat_ctx->Ku,-1,shell_mat_ctx->Mv); CHKERRQ(ierr); // f-Mrv
    ierr = PCApply(ctx->pC,shell_mat_ctx->Ku,shell_mat_ctx->u); CHKERRQ(ierr); //u = (aM+K)^(-1)(ru-Mrv)
    //VecView(shell_mat_ctx->u,PETSC_VIEWER_STDOUT_WORLD);
    //calculate velocities
    ierr = VecAXPY(shell_mat_ctx->v,shell_mat_ctx->ts_a,shell_mat_ctx->u); CHKERRQ(ierr); // v = v + a*u
    //VecView(shell_mat_ctx->v,PETSC_VIEWER_STDOUT_WORLD);
    //reverse
    ierr = VecZeroEntries(x); CHKERRQ(ierr);
    ierr = VecScatterBegin(shell_mat_ctx->scatterU,shell_mat_ctx->u,x,INSERT_VALUES,SCATTER_REVERSE); CHKERRQ(ierr);
    ierr = VecScatterEnd(shell_mat_ctx->scatterU,shell_mat_ctx->u,x,INSERT_VALUES,SCATTER_REVERSE); CHKERRQ(ierr);
    ierr = VecScatterBegin(shell_mat_ctx->scatterV,shell_mat_ctx->v,x,INSERT_VALUES,SCATTER_REVERSE); CHKERRQ(ierr);
    ierr = VecScatterEnd(shell_mat_ctx->scatterV,shell_mat_ctx->v,x,INSERT_VALUES,SCATTER_REVERSE); CHKERRQ(ierr);
    ierr = VecAssemblyBegin(x); CHKERRQ(ierr);
    ierr = VecAssemblyEnd(x); CHKERRQ(ierr);
    PetscFunctionReturn(0);
  }

  struct ShellResidualElement: public FEMethod {
    FieldInterface &mField;
    ShellResidualElement(FieldInterface &m_field): mField(m_field) {}

    //variables bellow need to be set by user
    MatShellCtx *shellMatCtx; 					///< pointer to shell matrix

    PetscErrorCode preProcess() {
      PetscFunctionBegin;

      PetscErrorCode ierr;
      if(ts_ctx != CTX_TSSETIFUNCTION) {
	SETERRQ(PETSC_COMM_SELF,MOFEM_NOT_IMPLEMENTED,"It is used to residual of velocities");
      }
      if(!shellMatCtx->iNitialized) {
	ierr = shellMatCtx->iNit(); CHKERRQ(ierr);
      }
      ierr = VecScatterBegin(shellMatCtx->scatterU,ts_u_t,shellMatCtx->u,INSERT_VALUES,SCATTER_FORWARD); CHKERRQ(ierr);      
      ierr = VecScatterEnd(shellMatCtx->scatterU,ts_u_t,shellMatCtx->u,INSERT_VALUES,SCATTER_FORWARD); CHKERRQ(ierr);      
      ierr = VecScatterBegin(shellMatCtx->scatterV,ts_u,shellMatCtx->v,INSERT_VALUES,SCATTER_FORWARD); CHKERRQ(ierr);
      ierr = VecScatterEnd(shellMatCtx->scatterV,ts_u,shellMatCtx->v,INSERT_VALUES,SCATTER_FORWARD); CHKERRQ(ierr);
      ierr = VecAXPY(shellMatCtx->v,-1,shellMatCtx->u); CHKERRQ(ierr);
      ierr = VecScatterBegin(shellMatCtx->scatterV,shellMatCtx->v,ts_F,ADD_VALUES,SCATTER_REVERSE); CHKERRQ(ierr);
      ierr = VecScatterEnd(shellMatCtx->scatterV,shellMatCtx->v,ts_F,ADD_VALUES,SCATTER_REVERSE); CHKERRQ(ierr);
      //VecView(shellMatCtx->v,PETSC_VIEWER_STDOUT_WORLD);

      PetscFunctionReturn(0);
    }

    PetscErrorCode postProcess() {
      PetscFunctionBegin;

      /*PetscErrorCode ierr;
      if(ts_ctx != CTX_TSSETIFUNCTION) {
	SETERRQ(PETSC_COMM_SELF,MOFEM_NOT_IMPLEMENTED,"It is used to residual of velocities");
      }
      if(!shellMatCtx->iNitialized) {
	ierr = shellMatCtx->iNit(); CHKERRQ(ierr);
      }
      ierr = VecScatterBegin(shellMatCtx->scatterU,ts_F,shellMatCtx->Ku,INSERT_VALUES,SCATTER_FORWARD); CHKERRQ(ierr);      
      ierr = VecScatterEnd(shellMatCtx->scatterU,ts_F,shellMatCtx->Ku,INSERT_VALUES,SCATTER_FORWARD); CHKERRQ(ierr);      
      ierr = VecScatterBegin(shellMatCtx->scatterV,ts_F,shellMatCtx->Mv,INSERT_VALUES,SCATTER_FORWARD); CHKERRQ(ierr);      
      ierr = VecScatterEnd(shellMatCtx->scatterV,ts_F,shellMatCtx->Mv,INSERT_VALUES,SCATTER_FORWARD); CHKERRQ(ierr);      

      double nrm2_ku,nrm2_mv;
      ierr = VecNorm(shellMatCtx->Ku,NORM_INFINITY,&nrm2_ku); CHKERRQ(ierr);
      ierr = VecNorm(shellMatCtx->Mv,NORM_INFINITY,&nrm2_mv); CHKERRQ(ierr);
      PetscPrintf(mField.get_comm(),"nrm2 U = %6.4e nrm2 V = %6.4e scale = %6.4e\n",nrm2_ku,nrm2_mv,nrm2_ku/fmax(nrm2_mv,nrm2_ku));
      //shellMatCtx->scale = nrm2_ku/fmax(nrm2_mv,nrm2_ku);
      //ierr = VecScale(shellMatCtx->Mv,shellMatCtx->scale); CHKERRQ(ierr);
      ierr = VecScatterBegin(shellMatCtx->scatterV,shellMatCtx->Mv,ts_F,INSERT_VALUES,SCATTER_REVERSE); CHKERRQ(ierr);      
      ierr = VecScatterEnd(shellMatCtx->scatterV,shellMatCtx->Mv,ts_F,INSERT_VALUES,SCATTER_REVERSE); CHKERRQ(ierr); */

      PetscFunctionReturn(0);

    }

  };

  #ifdef __DIRICHLETBC_HPP__

  /** \brief blocked element/problem
    *
    * Blocked element run loops for different problem than TS problem. It is
    * used to calculate matrices of shell matrix.
    *
    */
  struct ShellMatrixElement: public FEMethod {
    FieldInterface &mField;
    ShellMatrixElement(FieldInterface &m_field): mField(m_field) {}

    typedef pair<string,FEMethod*> LoopPairType;
    typedef vector<LoopPairType > LoopsToDoType;
    LoopsToDoType loopK; 	///< methods to calculate K shell matrix
    LoopsToDoType loopM; 	///< methods to calculate M shell matrix 
    LoopsToDoType loopAuxM; 	///< methods to calculate derivatives of inertia forces over displacements shell matrix 

    //variables bellow need to be set by user
    string problemName; 					///< name of shell problem
    MatShellCtx *shellMatCtx; 					///< pointer to shell matrix
    SpatialPositionsBCFEMethodPreAndPostProc *dirihletBcPtr; 	///< boundary conditions

    PetscErrorCode preProcess() {
      PetscFunctionBegin;
      PetscErrorCode ierr;

      if(ts_ctx != CTX_TSSETIJACOBIAN) {
	SETERRQ(PETSC_COMM_SELF,MOFEM_NOT_IMPLEMENTED,"It is used to calculate shell matrix only");
      }

      shellMatCtx->ts_a = ts_a;
      dirihletBcPtr->copy_ts(*((TSMethod*)this)); //copy context for TSMethod   

      dirihletBcPtr->dIag = 1;
      dirihletBcPtr->ts_B = shellMatCtx->K;
      ierr = MatZeroEntries(shellMatCtx->K); CHKERRQ(ierr);
      ierr = mField.problem_basic_method_preProcess(problemName,*dirihletBcPtr); CHKERRQ(ierr);
      LoopsToDoType::iterator itk = loopK.begin();
      for(;itk!=loopK.end();itk++) {
	itk->second->copy_ts(*((TSMethod*)this));
	itk->second->ts_B = shellMatCtx->K;
	ierr = mField.loop_finite_elements(problemName,itk->first,*itk->second); CHKERRQ(ierr);
      }
      LoopsToDoType::iterator itam = loopAuxM.begin();
      for(;itam!=loopAuxM.end();itam++) {
	itam->second->copy_ts(*((TSMethod*)this));
	itam->second->ts_B = shellMatCtx->K;
	ierr = mField.loop_finite_elements(problemName,itam->first,*itam->second); CHKERRQ(ierr);
      }
      ierr = mField.problem_basic_method_postProcess(problemName,*dirihletBcPtr); CHKERRQ(ierr);
      ierr = MatAssemblyBegin(shellMatCtx->K,MAT_FINAL_ASSEMBLY); CHKERRQ(ierr);
      ierr = MatAssemblyEnd(shellMatCtx->K,MAT_FINAL_ASSEMBLY); CHKERRQ(ierr);

      dirihletBcPtr->dIag = 0;
      dirihletBcPtr->ts_B = shellMatCtx->M;
      ierr = MatZeroEntries(shellMatCtx->M); CHKERRQ(ierr);
      //ierr = mField.problem_basic_method_preProcess(problemName,*dirihletBcPtr); CHKERRQ(ierr);
      LoopsToDoType::iterator itm = loopM.begin();
      for(;itm!=loopM.end();itm++) {
	itm->second->copy_ts(*((TSMethod*)this));
	itm->second->ts_B = shellMatCtx->M;
	ierr = mField.loop_finite_elements(problemName,itm->first,*itm->second); CHKERRQ(ierr);
      }
      ierr = mField.problem_basic_method_postProcess(problemName,*dirihletBcPtr); CHKERRQ(ierr);
      ierr = MatAssemblyBegin(shellMatCtx->M,MAT_FINAL_ASSEMBLY); CHKERRQ(ierr);
      ierr = MatAssemblyEnd(shellMatCtx->M,MAT_FINAL_ASSEMBLY); CHKERRQ(ierr);

      //barK
      ierr = MatZeroEntries(shellMatCtx->barK); CHKERRQ(ierr);
      ierr = MatCopy(shellMatCtx->K,shellMatCtx->barK,SAME_NONZERO_PATTERN); CHKERRQ(ierr);
      ierr = MatAXPY(shellMatCtx->barK,ts_a,shellMatCtx->M,SAME_NONZERO_PATTERN); CHKERRQ(ierr);
      ierr = MatAssemblyBegin(shellMatCtx->barK,MAT_FINAL_ASSEMBLY); CHKERRQ(ierr);
      ierr = MatAssemblyEnd(shellMatCtx->barK,MAT_FINAL_ASSEMBLY); CHKERRQ(ierr);

      //Matrix View
      //MatView(shellMatCtx->barK,PETSC_VIEWER_DRAW_WORLD);//PETSC_VIEWER_STDOUT_WORLD);
      //std::string wait;
      //std::cin >> wait;

      PetscFunctionReturn(0);
    } 

  };

  #endif //__DIRICHLETBC_HPP__

};


#endif //__CONVECTIVE_MASS_ELEMENT_HPP

/***************************************************************************//**
 * \defgroup convective_mass_elem Mass Element
 * \ingroup mofem_forces_and_sources 
 ******************************************************************************/



