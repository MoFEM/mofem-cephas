/* Copyright (C) 2013, Lukasz Kaczmarczyk (likask AT wp.pl)
 * --------------------------------------------------------------
 * DESCRIPTION:Analytical Dirichlet BC for Helmholtz operator and
 * The L^2 Norm of approximation error implemented.
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

#ifndef __ANALYTICALDIRIHLET_HELMHOLTZ_BC_HPP__
#define __ANALYTICALDIRIHLET_HELMHOLTZ_BC_HPP__

using namespace boost::numeric;
using namespace MoFEM;

struct AnalyticalDirihletBC {

  //ubals::vector<double> analyticalFunction(
    //double x,double y,double z) = 0;

  /** \brief finite element to appeximate analytical solution on surface
    */
  struct ApproxField {

    struct MyTriFE: public TriElementForcesAndSurcesCore {
      MyTriFE(FieldInterface &m_field): TriElementForcesAndSurcesCore(m_field) {}
      int getRule(int order) { return order; };
    };


    MyTriFE feApproxTri; 
    MyTriFE& getLoopFeApproxTri() { return feApproxTri; } 

	
	
	ApproxField(FieldInterface &m_field): feApproxTri(m_field) {};
	
    ublas::matrix<double> hoCoordsTri;

	//get the higher order approximation coordinates.
    struct OpHoCoordTri: public TriElementForcesAndSurcesCore::UserDataOperator {

      ublas::matrix<double> &hoCoordsTri;
      OpHoCoordTri(const string field_name,ublas::matrix<double> &ho_coords): 
	TriElementForcesAndSurcesCore::UserDataOperator(field_name),
	hoCoordsTri(ho_coords) {}

	/*	
	Cartesian coordinates for gaussian points inside elements
	X^coordinates = DOF dot* N
    */
      PetscErrorCode doWork(
	int side,EntityType type,DataForcesAndSurcesCore::EntData &data) {
	PetscFunctionBegin;

	try {

	  if(data.getFieldData().size()==0) PetscFunctionReturn(0);

	  hoCoordsTri.resize(data.getN().size1(),3);
	  if(type == MBVERTEX) {
	    hoCoordsTri.clear();
	  }

	  int nb_dofs = data.getFieldData().size();
	  //int nb_dofs = data.getN().size2();
	  for(unsigned int gg = 0;gg<data.getN().size1();gg++) {
		  
		  //cout << "\n nb_dofs = " << nb_dofs << "\n data.getN(gg) = " << data.getN(gg) << "\n data.getFieldData() =" << data.getFieldData() << std::endl;
	    for(int dd = 0;dd<3;dd++) {
		  hoCoordsTri(gg,dd) += cblas_ddot(nb_dofs/3,&data.getN(gg)[0],1,&data.getFieldData()[dd],3);
		  //std::cout <<  "\n data.getFieldData()[dd] =" << data.getFieldData()[dd] << std::endl;
	      //hoCoordsTri(gg,dd) += cblas_ddot(nb_dofs,&data.getN(gg)[0],1,&data.getFieldData()[dd],3); //calculate x,y,z in each GaussPts
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
 	
    /** \brief Lhs operaetar for triangle used to build matrix
      */
    struct OpLhsTri:public TriElementForcesAndSurcesCore::UserDataOperator {
		
		Mat B;
		bool solveBc;
		ublas::matrix<double> &hoCoords;
		OpLhsTri(const string field_name,ublas::matrix<double> &ho_coords,Mat _B): 
			TriElementForcesAndSurcesCore::UserDataOperator(field_name),
			hoCoords(ho_coords),B(_B),solveBc(false) { }
		OpLhsTri(const string field_name,ublas::matrix<double> &ho_coords): 
			TriElementForcesAndSurcesCore::UserDataOperator(field_name),
			hoCoords(ho_coords),solveBc(true) { }
		
		ublas::matrix<FieldData> NTN;
	  
	  /*	
	  Lhs mass matrix
	  A = N^T N
	  */
	  
      PetscErrorCode doWork(
	int row_side,int col_side,
	EntityType row_type,EntityType col_type,
	DataForcesAndSurcesCore::EntData &row_data,
	DataForcesAndSurcesCore::EntData &col_data) {
	PetscFunctionBegin;

	PetscErrorCode ierr;

	try {

	  if(row_data.getIndices().size()==0) PetscFunctionReturn(0);
	  if(col_data.getIndices().size()==0) PetscFunctionReturn(0);
  
	  unsigned int nb_row = row_data.getIndices().size();
	  unsigned int nb_col = col_data.getIndices().size();
	  //unsigned int nb_row = row_data.getN().size2();
	  //unsigned int nb_col = col_data.getN().size2();

	  if(nb_row != row_data.getIndices().size()) {
	    SETERRQ(PETSC_COMM_SELF,MOFEM_NOT_IMPLEMENTED,
	      "currently works only for scalar fields, extension to fields with higher rank need to be implemented");
	  }

	  NTN.resize(nb_row,nb_col);

	  for(unsigned int gg = 0;gg<row_data.getN().size1();gg++) {
	    double val = getGaussPts()(2,gg);
	    if(hoCoords.size1() == row_data.getN().size1()) {
	      double area = norm_2(getNormals_at_GaussPt(gg))*0.5; 
	      val *= area;
	    } else {
	      val *= getArea();
	    }
         
	    NTN.clear();
	    cblas_dger(CblasRowMajor,nb_row,nb_col,val,
	      &row_data.getN(gg)[0],1,&col_data.getN(gg)[0],1,
	      &NTN(0,0),nb_col);
		
        if(solveBc) {
	    ierr = MatSetValues(
	      (getFEMethod()->snes_B),
	      nb_row,&row_data.getIndices()[0],
	      nb_col,&col_data.getIndices()[0],
	      &NTN(0,0),ADD_VALUES); CHKERRQ(ierr);
		
		} else if(!solveBc){
		ierr = MatSetValues(
				   B,
				   nb_row,&row_data.getIndices()[0],
				   nb_col,&col_data.getIndices()[0],
				   &NTN(0,0),ADD_VALUES); CHKERRQ(ierr);
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
	
	
    /** \brief Rhs Triangle operaetar used to build matrix
      */
	
	struct OpRhsTri:public TriElementForcesAndSurcesCore::UserDataOperator {
		
		ublas::matrix<double> &hoCoords;
		double (*fUN)(double x,double y,double z,bool use_real);
		Vec C;
		bool solveBc;
		bool use_real;
		string fieldType;
		
		
		OpRhsTri(const string field_name,ublas::matrix<double> &ho_coords,
			  double (*fun)(double x,double y,double z,bool use_real),Vec _C): 
			TriElementForcesAndSurcesCore::UserDataOperator(field_name),
			hoCoords(ho_coords),fUN(fun),C(_C),fieldType(field_name),solveBc(false)  {}
		OpRhsTri(const string field_name,ublas::matrix<double> &ho_coords,
				 double (*fun)(double x,double y,double z,bool use_real)): 
			TriElementForcesAndSurcesCore::UserDataOperator(field_name),
			hoCoords(ho_coords),fUN(fun),fieldType(field_name),solveBc(true)  {}
		
		ublas::vector<FieldData> NTf;
	  
	  /*	
	  Rhs force vector merely with Dirichlet values
	  F = int_S N^T Fun dS
	  */
	  
      PetscErrorCode doWork(
	int side,EntityType type,DataForcesAndSurcesCore::EntData &data) {
	PetscFunctionBegin;
	PetscErrorCode ierr;
  
	try {

	  if(data.getIndices().size()==0) PetscFunctionReturn(0);
  
	  unsigned int nb_row = data.getIndices().size();
	  //unsigned int nb_row = data.getFieldData().size();
	  //unsigned int nb_row = data.getN().size2();
	  if(nb_row != data.getIndices().size()) {
	    SETERRQ(PETSC_COMM_SELF,MOFEM_NOT_IMPLEMENTED,
	      "currently works only for scalar fields, extension to fields with higher rank need to be implemented");
	  }
	  NTf.resize(nb_row);

	  for(unsigned int gg = 0;gg<data.getN().size1();gg++) {

	    double x,y,z;
	    double val = getGaussPts()(2,gg);
		
	    if(hoCoords.size1() == data.getN().size1()) {
	      double area = norm_2(getNormals_at_GaussPt(gg))*0.5; 
	      val *= area;
	      x = hoCoords(gg,0);
	      y = hoCoords(gg,1);
	      z = hoCoords(gg,2);
	    } else {
	      val *= getArea();
	      x = getCoordsAtGaussPts()(gg,0);
	      y = getCoordsAtGaussPts()(gg,1);
	      z = getCoordsAtGaussPts()(gg,2);
	    }
	    
		if(fieldType.compare(0,6,"rePRES")==0) {
			use_real = true;
		} else if(fieldType.compare(0,6,"imPRES")==0) {
			use_real = false;
		}
				
	    double a = fUN(x,y,z,use_real);
		
		
		
	    noalias(NTf) = data.getN(gg,nb_row)*a*val;
		
        if(solveBc){
	    ierr = VecSetValues(getFEMethod()->snes_f,data.getIndices().size(),
	      &data.getIndices()[0],&*NTf.data().begin(),ADD_VALUES); CHKERRQ(ierr);
		} else if(!solveBc){
		
	    ierr = VecSetValues(C,data.getIndices().size(),
							&data.getIndices()[0],&*NTf.data().begin(),ADD_VALUES); CHKERRQ(ierr);
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
	
  };

  struct DirichletBC: public DisplacementBCFEMethodPreAndPostProc {

    DirichletBC(
      FieldInterface& m_field,const string &field,Mat A,Vec X,Vec F): 
      DisplacementBCFEMethodPreAndPostProc(m_field,field,A,X,F),tRis_ptr(NULL) {}
    DirichletBC(
      FieldInterface& m_field,const string &field): 
      DisplacementBCFEMethodPreAndPostProc(m_field,field),tRis_ptr(NULL) {}

    Range *tRis_ptr;

    PetscErrorCode iNitalize() {
      PetscFunctionBegin;
      if(map_zero_rows.empty()) {
	if(tRis_ptr == NULL) {
	  SETERRQ(PETSC_COMM_SELF,1,"need to inicialsised from AnalyticalDirihletBC::solveProblem");
	}
	ierr = iNitalize(*tRis_ptr); CHKERRQ(ierr);
      }
      PetscFunctionReturn(0);
    }
  
    PetscErrorCode iNitalize(Range &tris) {
      PetscFunctionBegin;
      ParallelComm* pcomm = ParallelComm::get_pcomm(&mField.get_moab(),MYPCOMM_INDEX);
      Range ents;
      rval = mField.get_moab().get_connectivity(tris,ents,true); CHKERR_PETSC(rval);
      ierr = mField.get_moab().get_adjacencies(tris,1,false,ents,Interface::UNION); CHKERRQ(ierr); //2rd input mean the dim, 1 is edge.
      ents.merge(tris);
      for(Range::iterator eit = ents.begin();eit!=ents.end();eit++) {
    	for(_IT_NUMEREDDOFMOFEMENTITY_ROW_BY_NAME_ENT_PART_FOR_LOOP_(problemPtr,fieldName,*eit,pcomm->rank(),dof)) {
    	  map_zero_rows[dof->get_petsc_gloabl_dof_idx()] = dof->get_FieldData();
        }
      }
      dofsIndices.resize(map_zero_rows.size());
      dofsValues.resize(map_zero_rows.size());
      int ii = 0;
      map<DofIdx,FieldData>::iterator mit = map_zero_rows.begin();
      for(;mit!=map_zero_rows.end();mit++,ii++) { 
	dofsIndices[ii] = mit->first;
	dofsValues[ii] = mit->second;
      }
      PetscFunctionReturn(0);
    }


  };
  
  

  
  
  ApproxField approxField;
  Range tRis;
  AnalyticalDirihletBC(FieldInterface& m_field,Range &bc_tris): approxField(m_field),tRis(bc_tris) {};
  
  PetscErrorCode setApproxOps(
	  FieldInterface &m_field,
	  string re_field_name,
	  double (*fun)(double x,double y,double z,bool use_real),
	  string nodals_positions = "MESH_NODE_POSITIONS") {
	  PetscFunctionBegin;
	  if(m_field.check_field(nodals_positions)) {
		  approxField.getLoopFeApproxTri().get_op_to_do_Rhs().push_back(new ApproxField::OpHoCoordTri(nodals_positions,approxField.hoCoordsTri));
	  }
	  //loop over triangles
	  approxField.getLoopFeApproxTri().get_op_to_do_Rhs().push_back(new ApproxField::OpRhsTri(re_field_name,approxField.hoCoordsTri,fun));
	  approxField.getLoopFeApproxTri().get_op_to_do_Lhs().push_back(new ApproxField::OpLhsTri(re_field_name,approxField.hoCoordsTri));
	  //loop over tets
	  PetscFunctionReturn(0);
  }
  

  PetscErrorCode initializeBcProblem(
    FieldInterface &m_field,
    string problem,string fe,string re_field_name,
    string nodals_positions = "MESH_NODE_POSITIONS") {
    PetscFunctionBegin;
    PetscErrorCode ierr;
	//Add triangle elements
    ierr = m_field.add_finite_element(fe,MF_ZERO); CHKERRQ(ierr);
    ierr = m_field.modify_finite_element_add_field_row(fe,re_field_name); CHKERRQ(ierr);
    ierr = m_field.modify_finite_element_add_field_col(fe,re_field_name); CHKERRQ(ierr);
    ierr = m_field.modify_finite_element_add_field_data(fe,re_field_name); CHKERRQ(ierr);
    
	
    if(m_field.check_field(nodals_positions)) {
      ierr = m_field.modify_finite_element_add_field_data(fe,nodals_positions); CHKERRQ(ierr);
    }
    ierr = m_field.modify_problem_add_finite_element(problem,fe); CHKERRQ(ierr);
    ierr = m_field.add_ents_to_finite_element_by_TRIs(tRis,fe); CHKERRQ(ierr);
	
    PetscFunctionReturn(0);
  }
  

  Mat A;
  Vec D,F;
  KSP solver;
  
  PetscErrorCode setBcProblem(
    FieldInterface &m_field,string problem) {
    PetscFunctionBegin;
    PetscErrorCode ierr;

    ierr = m_field.VecCreateGhost(problem,ROW,&F); CHKERRQ(ierr);
    ierr = m_field.VecCreateGhost(problem,COL,&D); CHKERRQ(ierr);
    ierr = m_field.MatCreateMPIAIJWithArrays(problem,&A); CHKERRQ(ierr);

    ierr = KSPCreate(PETSC_COMM_WORLD,&solver); CHKERRQ(ierr);
    ierr = KSPSetOperators(solver,A,A); CHKERRQ(ierr);
    ierr = KSPSetFromOptions(solver); CHKERRQ(ierr);

    PetscFunctionReturn(0);
  }
  
  PetscErrorCode solveBcProblem(
    FieldInterface &m_field,string problem,string fe,DirichletBC &bc) {
    PetscFunctionBegin;
    PetscErrorCode ierr;

    ierr = VecZeroEntries(F); CHKERRQ(ierr);
    ierr = MatZeroEntries(A); CHKERRQ(ierr);

    approxField.getLoopFeApproxTri().snes_B = A;
    approxField.getLoopFeApproxTri().snes_f = F;
    ierr = m_field.loop_finite_elements(problem,fe,approxField.getLoopFeApproxTri()); CHKERRQ(ierr);

    ierr = VecGhostUpdateBegin(F,ADD_VALUES,SCATTER_REVERSE); CHKERRQ(ierr);
    ierr = VecGhostUpdateEnd(F,ADD_VALUES,SCATTER_REVERSE); CHKERRQ(ierr);
    ierr = VecAssemblyBegin(F); CHKERRQ(ierr);
    ierr = VecAssemblyEnd(F); CHKERRQ(ierr);
    ierr = MatAssemblyBegin(A,MAT_FINAL_ASSEMBLY); CHKERRQ(ierr);
    ierr = MatAssemblyEnd(A,MAT_FINAL_ASSEMBLY); CHKERRQ(ierr);
	
	int ii1,jj1;
	ierr=MatGetSize(A,&ii1,&jj1);
	
	//std::string wait;
	//ierr = MatView(A,PETSC_VIEWER_STDOUT_WORLD);
	//ierr = MatView(A,PETSC_VIEWER_DRAW_WORLD);
	//std::cin >> wait;
	
    ierr = KSPSolve(solver,F,D); CHKERRQ(ierr);
    ierr = VecGhostUpdateBegin(D,INSERT_VALUES,SCATTER_FORWARD); CHKERRQ(ierr);
    ierr = VecGhostUpdateEnd(D,INSERT_VALUES,SCATTER_FORWARD); CHKERRQ(ierr);
	
    ierr = m_field.set_global_VecCreateGhost(problem,ROW,D,INSERT_VALUES,SCATTER_REVERSE); CHKERRQ(ierr);

    bc.tRis_ptr = &tRis; 
    bc.map_zero_rows.clear();
    bc.dofsIndices.clear();
    bc.dofsValues.clear();
  
    PetscFunctionReturn(0);
  }

  PetscErrorCode destroyBcProblem() {
    PetscFunctionBegin;
    PetscErrorCode ierr;
    ierr = KSPDestroy(&solver); CHKERRQ(ierr);
    ierr = MatDestroy(&A); CHKERRQ(ierr);
    ierr = VecDestroy(&F); CHKERRQ(ierr);
    ierr = VecDestroy(&D); CHKERRQ(ierr);
    PetscFunctionReturn(0);
  }
  
};

#endif //__ANALYTICALDIRIHLET_HELMHOLTZ_BC_HPP__

