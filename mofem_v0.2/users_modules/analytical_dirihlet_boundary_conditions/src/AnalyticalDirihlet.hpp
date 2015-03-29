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

#ifndef __ANALYTICALDIRIHLETBC_HPP__
#define __ANALYTICALDIRIHLETBC_HPP__

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

    ApproxField(FieldInterface &m_field): feApprox(m_field) {};

    MyTriFE feApprox; 
    MyTriFE& getLoopFeApprox() { return feApprox; } 

    ublas::matrix<double> hoCoords;
    struct OpHoCoord: public TriElementForcesAndSurcesCore::UserDataOperator {

      ublas::matrix<double> &hoCoords;
      OpHoCoord(const string field_name,ublas::matrix<double> &ho_coords): 
	TriElementForcesAndSurcesCore::UserDataOperator(field_name),
	hoCoords(ho_coords) {}

      PetscErrorCode doWork(
	int side,EntityType type,DataForcesAndSurcesCore::EntData &data) {
	PetscFunctionBegin;

	try {

	  if(data.getFieldData().size()==0) PetscFunctionReturn(0);

	  hoCoords.resize(data.getN().size1(),3);
	  if(type == MBVERTEX) {
	    hoCoords.clear();
	  }

	  int nb_dofs = data.getFieldData().size();
	  for(unsigned int gg = 0;gg<data.getN().size1();gg++) {
	    for(int dd = 0;dd<3;dd++) {
	      hoCoords(gg,dd) += cblas_ddot(nb_dofs/3,&data.getN(gg)[0],1,&data.getFieldData()[dd],3);
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
  

    /** \brief Lhs operaetar used to build matrix
      */
    struct OpLhs:public TriElementForcesAndSurcesCore::UserDataOperator {

      ublas::matrix<double> &hoCoords;
      OpLhs(const string field_name,ublas::matrix<double> &ho_coords): 
	TriElementForcesAndSurcesCore::UserDataOperator(field_name),
	hoCoords(ho_coords) { }

      ublas::matrix<FieldData> NTN,transNTN;
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

	  if(nb_row != row_data.getIndices().size()) {
	    SETERRQ(PETSC_COMM_SELF,MOFEM_NOT_IMPLEMENTED,
	      "currently works only for scalar fields, extension to fields with higher rank need to be implemented");
	  }

	  NTN.resize(nb_row,nb_col);
	  NTN.clear();

	  for(unsigned int gg = 0;gg<row_data.getN().size1();gg++) {
	    double val = getGaussPts()(2,gg);
	    if(hoCoords.size1() == row_data.getN().size1()) {
        
	      double area = norm_2(getNormals_at_GaussPt(gg))*0.5; 

	      val *= area;
	    } else {
	      val *= getArea();
	    }


	    cblas_dger(CblasRowMajor,nb_row,nb_col,val,
	      &row_data.getN(gg)[0],1,&col_data.getN(gg)[0],1,
	      &NTN(0,0),nb_col);


	  }

	  ierr = MatSetValues(
	    (getFEMethod()->snes_B),
	    nb_row,&row_data.getIndices()[0],
	    nb_col,&col_data.getIndices()[0],
	    &NTN(0,0),ADD_VALUES); CHKERRQ(ierr);
	  if(row_side != col_side || row_type != col_type) {
	    transNTN.resize(nb_col,nb_row);
	    noalias(transNTN) = trans(NTN);
	    ierr = MatSetValues(
	      (getFEMethod()->snes_B),
	      nb_col,&col_data.getIndices()[0],
	      nb_row,&row_data.getIndices()[0],
	      &transNTN(0,0),ADD_VALUES); CHKERRQ(ierr);
	  }

	} catch (const std::exception& ex) {
	  ostringstream ss;
	  ss << "throw in method: " << ex.what() << endl;
	  SETERRQ(PETSC_COMM_SELF,1,ss.str().c_str());
	}

	PetscFunctionReturn(0);
      }
    };

    /** \brief Rhs operaetar used to build matrix
      */
    struct OpRhs:public TriElementForcesAndSurcesCore::UserDataOperator {

      ublas::matrix<double> &hoCoords;
      double (*fUN)(double x,double y,double z);
      OpRhs(const string field_name,ublas::matrix<double> &ho_coords,
	double (*fun)(double x,double y,double z)): 
	TriElementForcesAndSurcesCore::UserDataOperator(field_name),
	hoCoords(ho_coords),fUN(fun)  {}

      ublas::vector<FieldData> NTf;
      PetscErrorCode doWork(
	int side,EntityType type,DataForcesAndSurcesCore::EntData &data) {
	PetscFunctionBegin;
	PetscErrorCode ierr;
  
	try {

	  if(data.getIndices().size()==0) PetscFunctionReturn(0);
  
	  unsigned int nb_row = data.getIndices().size();
	  if(nb_row != data.getIndices().size()) {
	    SETERRQ(PETSC_COMM_SELF,MOFEM_NOT_IMPLEMENTED,
	      "currently works only for scalar fields, extension to fields with higher rank need to be implemented");
	  }
	  NTf.resize(nb_row);

	  for(unsigned int gg = 0;gg<data.getN().size1();gg++) {

	    double x,y,z;
	    double val = getGaussPts()(2,gg);
	    if(hoCoords.size1() == data.getN().size1()) {
	      double area = norm_2(getNormals_at_GaussPt(gg)); 
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
	    
	    double a = fUN(x,y,z);
	    noalias(NTf) = data.getN(gg,nb_row)*a*val;
	    ierr = VecSetValues(getFEMethod()->snes_f,data.getIndices().size(),
	      &data.getIndices()[0],&*NTf.data().begin(),ADD_VALUES); CHKERRQ(ierr);

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

  struct DirihletBC: public DisplacementBCFEMethodPreAndPostProc {

    DirihletBC(
      FieldInterface& m_field,const string &field,Mat A,Vec X,Vec F): 
      DisplacementBCFEMethodPreAndPostProc(m_field,field,A,X,F),tRis_ptr(NULL) {}
    DirihletBC(
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
      ierr = mField.get_moab().get_adjacencies(tris,1,false,ents,Interface::UNION); CHKERRQ(ierr);
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
    string field_name,
    double (*fun)(double x,double y,double z),
    string nodals_positions = "MESH_NODE_POSITIONS") {
    PetscFunctionBegin;
    if(m_field.check_field(nodals_positions)) {
	approxField.getLoopFeApprox().get_op_to_do_Rhs().push_back(new ApproxField::OpHoCoord(nodals_positions,approxField.hoCoords));
    }
    approxField.getLoopFeApprox().get_op_to_do_Rhs().push_back(new ApproxField::OpRhs(field_name,approxField.hoCoords,fun));
    approxField.getLoopFeApprox().get_op_to_do_Lhs().push_back(new ApproxField::OpLhs(field_name,approxField.hoCoords));
    PetscFunctionReturn(0);
  }

  PetscErrorCode initializeProblem(
    FieldInterface &m_field,
    string problem,string fe,string field,
    string nodals_positions = "MESH_NODE_POSITIONS") {
    PetscFunctionBegin;
    PetscErrorCode ierr;
    ierr = m_field.add_finite_element(fe,MF_ZERO); CHKERRQ(ierr);
    ierr = m_field.modify_finite_element_add_field_row(fe,field); CHKERRQ(ierr);
    ierr = m_field.modify_finite_element_add_field_col(fe,field); CHKERRQ(ierr);
    ierr = m_field.modify_finite_element_add_field_data(fe,field); CHKERRQ(ierr);
    ierr = m_field.modify_finite_element_add_field_data(fe,field); CHKERRQ(ierr);
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
  PetscErrorCode setProblem(
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
  
  PetscErrorCode solveProblem(
    FieldInterface &m_field,string problem,string fe,DirihletBC &bc) {
    PetscFunctionBegin;
    PetscErrorCode ierr;

    ierr = VecZeroEntries(F); CHKERRQ(ierr);
    ierr = MatZeroEntries(A); CHKERRQ(ierr);

    approxField.getLoopFeApprox().snes_B = A;
    approxField.getLoopFeApprox().snes_f = F;
    ierr = m_field.loop_finite_elements(problem,fe,approxField.getLoopFeApprox()); CHKERRQ(ierr);
    ierr = VecGhostUpdateBegin(F,ADD_VALUES,SCATTER_REVERSE); CHKERRQ(ierr);
    ierr = VecGhostUpdateEnd(F,ADD_VALUES,SCATTER_REVERSE); CHKERRQ(ierr);
    ierr = VecAssemblyBegin(F); CHKERRQ(ierr);
    ierr = VecAssemblyEnd(F); CHKERRQ(ierr);
    ierr = MatAssemblyBegin(A,MAT_FINAL_ASSEMBLY); CHKERRQ(ierr);
    ierr = MatAssemblyEnd(A,MAT_FINAL_ASSEMBLY); CHKERRQ(ierr);

    ierr = KSPSolve(solver,F,D); CHKERRQ(ierr);
    ierr = VecGhostUpdateBegin(D,INSERT_VALUES,SCATTER_FORWARD); CHKERRQ(ierr);
    ierr = VecGhostUpdateEnd(D,INSERT_VALUES,SCATTER_FORWARD); CHKERRQ(ierr);

    ierr = m_field.set_global_ghost_vector(problem,ROW,D,INSERT_VALUES,SCATTER_REVERSE); CHKERRQ(ierr);

    bc.tRis_ptr = &tRis; 
    bc.map_zero_rows.clear();
    bc.dofsIndices.clear();
    bc.dofsValues.clear();
  
    PetscFunctionReturn(0);
  }

  PetscErrorCode destroyProblem() {
    PetscFunctionBegin;
    PetscErrorCode ierr;
    ierr = KSPDestroy(&solver); CHKERRQ(ierr);
    ierr = MatDestroy(&A); CHKERRQ(ierr);
    ierr = VecDestroy(&F); CHKERRQ(ierr);
    ierr = VecDestroy(&D); CHKERRQ(ierr);
    PetscFunctionReturn(0);
  }



};

#endif //__ANALYTICALDIRIHLETBC_HPP__

