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

struct CohesiveInterfaceElement {

  struct CommonData {
    ublas::matrix<double> gapGlob;
    ublas::matrix<double> gapLoc;
    ublas::vector<ublas::matrix<double> > R;
  };
  CommonData commonData;

  struct MyPrism: public MoFEM::FlatPrismElementForcesAndSurcesCore {
    MyPrism(FieldInterface &m_field): MoFEM::FlatPrismElementForcesAndSurcesCore(m_field) {}
    int getRule(int order) { return 10; };
  };
  MyPrism feRhs;
  MyPrism feLhs;
  MyPrism feHistory;

  CohesiveInterfaceElement(FieldInterface &m_field):
    feRhs(m_field),feLhs(m_field),feHistory(m_field) {};

  MyPrism& getFeRhs() { return feRhs; }
  MyPrism& getFeLhs() { return feLhs; }
  MyPrism& getFeHistory() { return feHistory; }

  /** \brief constritutive (physical) equation for interface
    */
  struct PhysicalEquation {

    FieldInterface &mField;
    bool isInitialised;
    PhysicalEquation(FieldInterface &m_field): 
      mField(m_field),isInitialised(false) {};

    double h,youngModulus,beta,ft,Gf;
    Range pRisms;
    Tag thKappa,thDamagedPrism;

    double E0,g0,kappa1;
    PetscErrorCode iNitailise(const FEMethod *fe_method) {
      PetscFunctionBegin;
      ErrorCode rval;
      double def_damaged = 0;
      rval = mField.get_moab().tag_get_handle(
	"DAMAGED_PRISM",1,MB_TYPE_INTEGER,thDamagedPrism,MB_TAG_CREAT|MB_TAG_SPARSE,&def_damaged); CHKERR_THROW(rval);
      const int def_len = 0;
      rval = mField.get_moab().tag_get_handle("_KAPPA",def_len,MB_TYPE_DOUBLE,
	thKappa,MB_TAG_CREAT|MB_TAG_SPARSE|MB_TAG_VARLEN,NULL); CHKERR_PETSC(rval);
      E0 = youngModulus/h;
      g0 = ft/E0;
      kappa1 = 2*Gf/ft;
      PetscFunctionReturn(0);
    }

    double calcG(int gg,ublas::matrix<double> gap_loc) {
      return sqrt(pow(gap_loc(gg,0),2)+beta*(pow(gap_loc(gg,1),2)+pow(gap_loc(gg,2),2)));
    }

    double *kappaPtr;
    int kappaSize;
    PetscErrorCode getKappa(int nb_gauss_pts,const FEMethod *fe_method) {
      PetscFunctionBegin;
      EntityHandle ent = fe_method->fePtr->get_ent();
      ErrorCode rval;
      rval = mField.get_moab().tag_get_by_ptr(thKappa,&ent,1,(const void **)&kappaPtr,&kappaSize); 
      if(rval != MB_SUCCESS || kappaSize != nb_gauss_pts) {
	ublas::vector<double> kappa;
	kappa.resize(nb_gauss_pts);
	kappa.clear();
	int tag_size[1]; 
	tag_size[0] = nb_gauss_pts;
	void const* tag_data[] = { &kappa[0] };
	rval = mField.get_moab().tag_set_by_ptr(thKappa,&ent,1,tag_data,tag_size);  CHKERR_PETSC(rval);
	rval = mField.get_moab().tag_get_by_ptr(thKappa,&ent,1,(const void **)&kappaPtr,&kappaSize);  CHKERR_PETSC(rval);
      }
      PetscFunctionReturn(0);
    }

    ublas::matrix<double> Dglob,Dloc;
    PetscErrorCode calcDglob(const double omega,ublas::matrix<double> &R) {
      PetscFunctionBegin;
      Dglob.resize(3,3);
      Dloc.resize(3,3);
      Dloc.clear();
      double E = (1-omega)*E0;
      Dloc(0,0) = E;
      Dloc(1,1) = E;
      Dloc(2,2) = E;
      Dglob = prod( Dloc, R );
      Dglob = prod( trans(R), Dglob );
      PetscFunctionReturn(0);
    }

    PetscErrorCode calcOmega(const double kappa,double& omega) {
      PetscFunctionBegin;
      omega = 0;
      if(kappa>=kappa1) {
	omega = 1;
	PetscFunctionReturn(0);
      } else if(kappa>0) {
	double a = (2.0*Gf*E0+ft*ft)*kappa;
	double b = (ft+E0*kappa)*Gf;
	omega = 0.5*a/b;
      }
      PetscFunctionReturn(0);
    }

    PetscErrorCode calcTangetDglob(const double omega,double g,const ublas::vector<double>& gap_loc,ublas::matrix<double> &R) {
      PetscFunctionBegin;
      Dglob.resize(3,3);
      Dloc.resize(3,3);
      double domega = 0.5*(2*Gf*E0+ft*ft)/((ft+(g-ft/E0)*E0)*Gf) - 0.5*((g-ft/E0)*(2*Gf*E0+ft*ft)*E0)/(pow(ft+(g-ft/E0)*E0,2)*Gf);
      Dloc.resize(3,3);
      //r0
      Dloc(0,0) = (1-omega)*E0 - domega*E0*gap_loc[0]*gap_loc[0]/g;
      Dloc(0,1) = -domega*E0*gap_loc[0]*beta*gap_loc[1]/g;
      Dloc(0,2) = -domega*E0*gap_loc[0]*beta*gap_loc[2]/g;
      //r1
      Dloc(1,0) = -domega*E0*gap_loc[1]*gap_loc[0]/g;
      Dloc(1,1) = (1-omega)*E0 - domega*E0*gap_loc[1]*beta*gap_loc[1]/g;
      Dloc(1,2) = -domega*E0*gap_loc[1]*beta*gap_loc[2]/g;
      //r2
      Dloc(2,0) = -domega*E0*gap_loc[2]*gap_loc[0]/g;
      Dloc(2,1) = -domega*E0*gap_loc[2]*beta*gap_loc[1]/g;
      Dloc(2,2) = (1-omega)*E0 - domega*E0*gap_loc[2]*beta*gap_loc[2]/g;
      Dglob = prod(Dloc,R);
      Dglob = prod(trans(R),Dglob);
      PetscFunctionReturn(0);
    }
     
    virtual PetscErrorCode calculateTraction(
      ublas::vector<double> &traction,
      int gg,CommonData &common_data,
      const FEMethod *fe_method) {
      PetscFunctionBegin;
      PetscErrorCode ierr;
      if(!isInitialised) {
	ierr = iNitailise(fe_method); CHKERRQ(ierr);
	isInitialised = true;
      }
      if(gg==0) {
	ierr = getKappa(common_data.gapGlob.size1(),fe_method); CHKERRQ(ierr);
      }
      double g = calcG(gg,common_data.gapLoc);
      double kappa = fmax(g-g0,kappaPtr[gg]);
      double omega = 0;
      ierr = calcOmega(kappa,omega); CHKERRQ(ierr);
      //cerr << gg << " " << omega << endl;
      ierr = calcDglob(omega,common_data.R[gg]); CHKERRQ(ierr);
      traction.resize(3);
      ublas::matrix_row<ublas::matrix<double> > gap_glob(common_data.gapGlob,gg);
      noalias(traction) = prod(Dglob,gap_glob);
      PetscFunctionReturn(0);
    }

    virtual PetscErrorCode calculateTangentStiffeness(
      ublas::matrix<double> &tangent_matrix,
      int gg,CommonData &common_data,
      const FEMethod *fe_method) {
      PetscFunctionBegin;
      PetscErrorCode ierr;
      try {
      if(!isInitialised) {
	ierr = iNitailise(fe_method); CHKERRQ(ierr);
	isInitialised = true;
      }
      if(gg==0) {
	ierr = getKappa(common_data.gapGlob.size1(),fe_method); CHKERRQ(ierr);
      }
      double g = calcG(gg,common_data.gapLoc);
      double kappa = fmax(g-g0,kappaPtr[gg]);
      double omega = 0;
      ierr = calcOmega(kappa,omega); CHKERRQ(ierr);
      //cerr << gg << " " << omega << endl;
      int iter;
      ierr = SNESGetIterationNumber(fe_method->snes,&iter); CHKERRQ(ierr);
      if((kappa <= kappaPtr[gg])||(kappa>=kappa1)||(iter <= 1)) {
	ierr = calcDglob(omega,common_data.R[gg]); CHKERRQ(ierr);
      } else {
	ublas::matrix_row<ublas::matrix<double> > g_loc(common_data.gapLoc,gg);
	ierr = calcTangetDglob(omega,g,g_loc,common_data.R[gg]);
      }
      tangent_matrix.resize(3,3);
      noalias(tangent_matrix) = Dglob;
      //cerr << "t " << tangent_matrix << endl;
      } catch (const std::exception& ex) {
	ostringstream ss;
	ss << "throw in method: " << ex.what() << endl;
	SETERRQ(PETSC_COMM_SELF,1,ss.str().c_str());
      }
      PetscFunctionReturn(0);
    }

    virtual PetscErrorCode updateHistory(
      CommonData &common_data,const FEMethod *fe_method) {
      PetscFunctionBegin;
      ErrorCode rval;
      PetscErrorCode ierr;
      if(!isInitialised) {
	ierr = iNitailise(fe_method); CHKERRQ(ierr);
	isInitialised = true;
      }
      ierr = getKappa(common_data.gapGlob.size1(),fe_method); CHKERRQ(ierr);
      bool all_gauss_pts_damaged = true;
      for(unsigned int gg = 0;gg<common_data.gapGlob.size1();gg++) {
	double omega = 0;
	double g = calcG(gg,common_data.gapLoc);
	double kappa = fmax(g-g0,kappaPtr[gg]);
	kappaPtr[gg] = kappa;
	ierr = calcOmega(kappa,omega); CHKERRQ(ierr);
	//if(omega < 1.) {
	  all_gauss_pts_damaged = false;
	//}
      }
      if(all_gauss_pts_damaged) {
	EntityHandle ent = fe_method->fePtr->get_ent();
	int set_prism_as_demaged = 1;
	rval = mField.get_moab().tag_set_data(thDamagedPrism,&ent,1,&set_prism_as_demaged); CHKERR_PETSC(rval);
      }
      PetscFunctionReturn(0);
    }

  };
  

  /** \brief Set negative sign to shape functions on face 4
    */
  struct OpSetSignToShapeFunctions: public FlatPrismElementForcesAndSurcesCore::UserDataOperator {

    OpSetSignToShapeFunctions(const string field_name):
      FlatPrismElementForcesAndSurcesCore::UserDataOperator(field_name) {}

    PetscErrorCode doWork(
      int side,EntityType type,DataForcesAndSurcesCore::EntData &data) {
      PetscFunctionBegin;
      if(data.getN().size1()==0)  PetscFunctionReturn(0);
      if(data.getN().size2()==0)  PetscFunctionReturn(0);
      switch(type) {
	case MBVERTEX:
	  for(unsigned int gg = 0;gg<data.getN().size1();gg++) {
	    for(int nn = 3;nn<6;nn++) {
	      data.getN()(gg,nn) *= -1;
	    }
	  }
	  break;
	case MBEDGE:
	  if(side < 3) PetscFunctionReturn(0);
	  data.getN() *= -1;
	  break;
	case MBTRI:
	  if(side == 3) PetscFunctionReturn(0);
	  data.getN() *= -1;
	  break;
	default:
	  SETERRQ(PETSC_COMM_SELF,1,"data inconsitency");
      }
      PetscFunctionReturn(0);
    }

  };

  struct OpCalculateGapGlobal: public FlatPrismElementForcesAndSurcesCore::UserDataOperator {

    CommonData &commonData;
    OpCalculateGapGlobal(const string field_name,CommonData &common_data):
      FlatPrismElementForcesAndSurcesCore::UserDataOperator(field_name),
      commonData(common_data) {}

    PetscErrorCode doWork(
      int side,EntityType type,DataForcesAndSurcesCore::EntData &data) {
      PetscFunctionBegin;
      try {
	int nb_dofs = data.getIndices().size();
	if(nb_dofs == 0) PetscFunctionReturn(0); 
	int nb_gauss_pts = data.getN().size1();
	if(type == MBVERTEX) {
	  commonData.R.resize(nb_gauss_pts);
	  for(int gg = 0;gg<nb_gauss_pts;gg++) {
	    commonData.R[gg].resize(3,3);
	    double nrm2_normal = 0;
	    double nrm2_tangent1 = 0;
	    double nrm2_tangent2 = 0;
	    for(int dd = 0;dd<3;dd++) {
	      nrm2_normal += pow(getNormals_at_GaussPtF3()(gg,dd),2);	
	      nrm2_tangent1 += pow(getTangent1_at_GaussPtF3()(gg,dd),2);
	      nrm2_tangent2 += pow(getTangent2_at_GaussPtF3()(gg,dd),2);
	    }
	    nrm2_normal = sqrt(nrm2_normal);
	    nrm2_tangent1 = sqrt(nrm2_tangent1);
	    nrm2_tangent2 = sqrt(nrm2_tangent2);
	    for(int dd = 0;dd<3;dd++) {
	      commonData.R[gg](0,dd) = getNormals_at_GaussPtF3()(gg,dd)/nrm2_normal;
	      commonData.R[gg](1,dd) = getTangent1_at_GaussPtF3()(gg,dd)/nrm2_tangent1;
	      commonData.R[gg](2,dd) = getTangent2_at_GaussPtF3()(gg,dd)/nrm2_tangent2;
	    }
	  }
	}
	if(type == MBVERTEX) {
	  commonData.gapGlob.resize(nb_gauss_pts,3);
	  commonData.gapGlob.clear();
	}
	for(int gg = 0;gg<nb_gauss_pts;gg++) {
	  for(int dd = 0;dd<3;dd++) {
	    commonData.gapGlob(gg,dd) += cblas_ddot(
	      nb_dofs/3,&data.getN(gg)[0],1,&data.getFieldData()[dd],3);
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

  struct OpCalculateGapLocal: public FlatPrismElementForcesAndSurcesCore::UserDataOperator {

    CommonData &commonData;
    OpCalculateGapLocal(const string field_name,CommonData &common_data):
      FlatPrismElementForcesAndSurcesCore::UserDataOperator(field_name),
      commonData(common_data) {}

    PetscErrorCode doWork(
      int side,EntityType type,DataForcesAndSurcesCore::EntData &data) {
      PetscFunctionBegin;
      try {
	if(type == MBVERTEX) {
	  int nb_gauss_pts = data.getN().size1();
	  commonData.gapLoc.resize(nb_gauss_pts,3);
	  for(int gg = 0;gg<nb_gauss_pts;gg++) {
	    ublas::matrix_row<ublas::matrix<double> > gap_glob(commonData.gapGlob,gg);
	    ublas::matrix_row<ublas::matrix<double> > gap_loc(commonData.gapLoc,gg);
	    gap_loc = prod(commonData.R[gg],gap_glob);
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

  struct OpRhs: public FlatPrismElementForcesAndSurcesCore::UserDataOperator {

    CommonData &commonData;
    PhysicalEquation &physicalEqations;
    OpRhs(const string field_name,CommonData &common_data,PhysicalEquation &physical_eqations):
      FlatPrismElementForcesAndSurcesCore::UserDataOperator(field_name),
      commonData(common_data),physicalEqations(physical_eqations) {}

    ublas::vector<double> traction,Nf;
    PetscErrorCode doWork(
      int side,EntityType type,DataForcesAndSurcesCore::EntData &data) {
      PetscFunctionBegin;
      PetscErrorCode ierr;
      try {
	int nb_dofs = data.getIndices().size();
	if(nb_dofs == 0) PetscFunctionReturn(0); 
	if(physicalEqations.pRisms.find(getMoFEMFEPtr()->get_ent()) == physicalEqations.pRisms.end()) {
	  PetscFunctionReturn(0);
	}
	Nf.resize(nb_dofs);
	Nf.clear();
	int nb_gauss_pts = data.getN().size1();
	for(int gg = 0;gg<nb_gauss_pts;gg++) {
	  ierr = physicalEqations.calculateTraction(traction,gg,commonData,getFEMethod()); CHKERRQ(ierr);
	  double w = getGaussPts()(2,gg)*cblas_dnrm2(3,&getNormals_at_GaussPtF3()(gg,0),1)*0.5;
	  for(int nn = 0;nn<nb_dofs/3;nn++) {
	    for(int dd = 0;dd<3;dd++) {
	      Nf[3*nn+dd] += w*data.getN(gg)[nn]*traction[dd];
	    }
	  }
	}
	ierr = VecSetValues(getFEMethod()->snes_f,
	  data.getIndices().size(),&data.getIndices()[0],&Nf[0],ADD_VALUES); CHKERRQ(ierr);
      } catch (const std::exception& ex) {
	ostringstream ss;
	ss << "throw in method: " << ex.what() << endl;
	SETERRQ(PETSC_COMM_SELF,1,ss.str().c_str());
      }
      PetscFunctionReturn(0);
    }

  };

  struct OpLhs: public FlatPrismElementForcesAndSurcesCore::UserDataOperator {

    CommonData &commonData;
    PhysicalEquation &physicalEqations;
    OpLhs(const string field_name,CommonData &common_data,PhysicalEquation &physical_eqations):
      FlatPrismElementForcesAndSurcesCore::UserDataOperator(field_name),
      commonData(common_data),physicalEqations(physical_eqations) { symm = false; }

    ublas::matrix<double> K,D,ND;
    PetscErrorCode doWork(
      int row_side,int col_side,
      EntityType row_type,EntityType col_type,
      DataForcesAndSurcesCore::EntData &row_data,
      DataForcesAndSurcesCore::EntData &col_data) {
      PetscFunctionBegin;
      PetscErrorCode ierr;
      try {
	int nb_row = row_data.getIndices().size();
	if(nb_row == 0) PetscFunctionReturn(0);
	int nb_col = col_data.getIndices().size();
	if(nb_col == 0) PetscFunctionReturn(0);
	if(physicalEqations.pRisms.find(getMoFEMFEPtr()->get_ent()) 
	  == physicalEqations.pRisms.end()) {
	  PetscFunctionReturn(0);
	}
	//cerr << row_side << " " << row_type << " " << row_data.getN() << endl;
	//cerr << col_side << " " << col_type << " " << col_data.getN() << endl;
	ND.resize(nb_row,3);
	K.resize(nb_row,nb_col);
	K.clear();
	int nb_gauss_pts = row_data.getN().size1();
	for(int gg = 0;gg<nb_gauss_pts;gg++) {
	  ierr = physicalEqations.calculateTangentStiffeness(D,gg,commonData,getFEMethod()); CHKERRQ(ierr);
	  double w = getGaussPts()(2,gg)*cblas_dnrm2(3,&getNormals_at_GaussPtF3()(gg,0),1)*0.5;
	  ND.clear();
	  for(int nn = 0; nn<nb_row/3;nn++) {
	    for(int dd = 0;dd<3;dd++) {
	      for(int DD = 0;DD<3;DD++) {
		ND(3*nn+dd,DD) += row_data.getN(gg)[nn]*D(dd,DD);
	      }
	    }
	  }
	  for(int nn = 0; nn<nb_row/3; nn++) {
	    for(int dd = 0;dd<3;dd++) {  
	      for(int NN = 0; NN<nb_col/3; NN++) {
		for(int DD = 0; DD<3;DD++) {
		  K(3*nn+dd,3*NN+DD) += w*ND(3*nn+dd,DD)*col_data.getN(gg)[NN];
		}
	      }
	    }
	  }
	}
	ierr = MatSetValues(getFEMethod()->snes_B,
	  nb_row,&row_data.getIndices()[0],
	  nb_col,&col_data.getIndices()[0],
	  &K(0,0),ADD_VALUES); CHKERRQ(ierr);
      } catch (const std::exception& ex) {
	ostringstream ss;
	ss << "throw in method: " << ex.what() << endl;
	SETERRQ(PETSC_COMM_SELF,1,ss.str().c_str());
      }
      PetscFunctionReturn(0);
    }

  };

  struct OpHistory: public FlatPrismElementForcesAndSurcesCore::UserDataOperator {

    CommonData &commonData;
    PhysicalEquation &physicalEqations;
    OpHistory(const string field_name,CommonData &common_data,PhysicalEquation &physical_eqations):
      FlatPrismElementForcesAndSurcesCore::UserDataOperator(field_name),
      commonData(common_data),physicalEqations(physical_eqations) {}

    PetscErrorCode doWork(
      int side,EntityType type,DataForcesAndSurcesCore::EntData &data) {
      PetscFunctionBegin;
      PetscErrorCode ierr;
      if(type != MBVERTEX) PetscFunctionReturn(0);
      if(physicalEqations.pRisms.find(getMoFEMFEPtr()->get_ent()) == physicalEqations.pRisms.end()) {
	PetscFunctionReturn(0);
      }
      ierr = physicalEqations.updateHistory(commonData,getFEMethod()); CHKERRQ(ierr);
      PetscFunctionReturn(0);
    }

  };

  PetscErrorCode addOps(const string field_name,boost::ptr_vector<CohesiveInterfaceElement::PhysicalEquation> &interfaces) {
    PetscFunctionBegin;

    //Rhs
    feRhs.get_op_to_do_Rhs().push_back(new OpSetSignToShapeFunctions(field_name));
    feRhs.get_op_to_do_Rhs().push_back(new OpCalculateGapGlobal(field_name,commonData));
    feRhs.get_op_to_do_Rhs().push_back(new OpCalculateGapLocal(field_name,commonData));
    //Lhs
    feLhs.get_op_to_do_Rhs().push_back(new OpSetSignToShapeFunctions(field_name));
    feLhs.get_op_to_do_Rhs().push_back(new OpCalculateGapGlobal(field_name,commonData));
    feLhs.get_op_to_do_Rhs().push_back(new OpCalculateGapLocal(field_name,commonData));
    //History
    feHistory.get_op_to_do_Rhs().push_back(new OpSetSignToShapeFunctions(field_name));
    feHistory.get_op_to_do_Rhs().push_back(new OpCalculateGapGlobal(field_name,commonData));
    feHistory.get_op_to_do_Rhs().push_back(new OpCalculateGapLocal(field_name,commonData));

    //add equations/data for physical inerfaces
    boost::ptr_vector<CohesiveInterfaceElement::PhysicalEquation>::iterator pit;
    for(pit = interfaces.begin();pit!=interfaces.end();pit++) {
      feRhs.get_op_to_do_Rhs().push_back(new OpRhs(field_name,commonData,*pit));
      feLhs.get_op_to_do_Lhs().push_back(new OpLhs(field_name,commonData,*pit));
      feHistory.get_op_to_do_Rhs().push_back(new OpHistory(field_name,commonData,*pit));
    }

    PetscFunctionReturn(0);
  }

};
