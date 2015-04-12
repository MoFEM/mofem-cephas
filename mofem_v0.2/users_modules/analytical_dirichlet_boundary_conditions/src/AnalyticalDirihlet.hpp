/** \file AnalyticalDirihlet.hpp

  Enforce Dirichlet boundary condition for given analytical function,

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

/** \brief Analytical Dirichlet boundary conditions
  \ingroup mofem_forces_and_sources
  */
struct AnalyticalDirihletBC {

  /** \brief finite element to appeximate analytical solution on surface
    */
  struct ApproxField {


    struct MyTriFE: public TriElementForcesAndSurcesCore {

      int addToRule; ///< this is add to integration rule if 2nd order geometry approximation
      MyTriFE(FieldInterface &m_field): TriElementForcesAndSurcesCore(m_field),addToRule(1) {}
      int getRule(int order) { return order+addToRule; };

    };

    ApproxField(FieldInterface &m_field): feApprox(m_field) {}
    virtual ~ApproxField() {}

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

      ublas::matrix<FieldData> NN,transNN;
      PetscErrorCode doWork(
	int row_side,int col_side,
	EntityType row_type,EntityType col_type,
	DataForcesAndSurcesCore::EntData &row_data,
	DataForcesAndSurcesCore::EntData &col_data) {
	PetscFunctionBegin;

	PetscFunctionBegin;

	if(row_data.getIndices().size()==0) PetscFunctionReturn(0);
	if(col_data.getIndices().size()==0) PetscFunctionReturn(0);

	PetscErrorCode ierr;

	const FENumeredDofMoFEMEntity *dof_ptr;
	ierr = getMoFEMFEPtr()->get_row_dofs_by_petsc_gloabl_dof_idx(row_data.getIndices()[0],&dof_ptr); CHKERRQ(ierr);
	int rank = dof_ptr->get_max_rank();

	int nb_row_dofs = row_data.getIndices().size()/rank;
	int nb_col_dofs = col_data.getIndices().size()/rank;

	NN.resize(nb_row_dofs,nb_col_dofs);
	NN.clear();
	  
	unsigned int nb_gauss_pts = row_data.getN().size1();
	for(unsigned int gg = 0;gg<nb_gauss_pts;gg++) {

	  double w = getGaussPts()(2,gg);
	  if(hoCoords.size1() == row_data.getN().size1()) {
	  
	    // higher order element
	    double area = norm_2(getNormals_at_GaussPt(gg))*0.5; 
	    w *= area;

	  } else {
	    
	    //linear element
	    w *= getArea();

	  }

	  cblas_dger(CblasRowMajor,
	    nb_row_dofs,nb_col_dofs,
	    w,&row_data.getN()(gg,0),1,&col_data.getN()(gg,0),1,
	    &*NN.data().begin(),nb_col_dofs);

	}
      
	if( (row_type != col_type) || (row_side != col_side) ) {
	  transNN.resize(nb_col_dofs,nb_row_dofs);
	  ublas::noalias(transNN) = trans(NN);
	}

	double *data = &*NN.data().begin();
	double *trans_data = &*transNN.data().begin();
	ublas::vector<DofIdx> row_indices,col_indices;
	row_indices.resize(nb_row_dofs);
	col_indices.resize(nb_col_dofs);

	for(int rr = 0;rr < rank; rr++) {
      
	  if((row_data.getIndices().size()%rank)!=0) {
	    SETERRQ(PETSC_COMM_SELF,MOFEM_DATA_INCONSISTENCT,"data inconsistency");
	  }

	  if((col_data.getIndices().size()%rank)!=0) {
	    SETERRQ(PETSC_COMM_SELF,MOFEM_DATA_INCONSISTENCT,"data inconsistency");
	  }

	  unsigned int nb_rows;
	  unsigned int nb_cols;
	  int *rows;
	  int *cols;


	  if(rank > 1) {

	    ublas::noalias(row_indices) = ublas::vector_slice<ublas::vector<DofIdx> >
	      (row_data.getIndices(), ublas::slice(rr, rank, row_data.getIndices().size()/rank));
	    ublas::noalias(col_indices) = ublas::vector_slice<ublas::vector<DofIdx> >
	      (col_data.getIndices(), ublas::slice(rr, rank, col_data.getIndices().size()/rank));

	    nb_rows = row_indices.size();
	    nb_cols = col_indices.size();
	    rows = &*row_indices.data().begin();
	    cols = &*col_indices.data().begin();

	  } else {

	    nb_rows = row_data.getIndices().size();
	    nb_cols = col_data.getIndices().size();
	    rows = &*row_data.getIndices().data().begin();
	    cols = &*col_data.getIndices().data().begin();

	  }

	  if(nb_rows != NN.size1()) {
	    SETERRQ(PETSC_COMM_SELF,MOFEM_DATA_INCONSISTENCT,"data inconsistency");
	  } 
	  if(nb_cols != NN.size2()) {
	    SETERRQ(PETSC_COMM_SELF,MOFEM_DATA_INCONSISTENCT,"data inconsistency");
	  } 
  
	  ierr = MatSetValues(getFEMethod()->snes_B,nb_rows,rows,nb_cols,cols,data,ADD_VALUES); CHKERRQ(ierr);
	  if( (row_type != col_type) || (row_side != col_side) ) {
	    if(nb_rows != transNN.size2()) {
	      SETERRQ(PETSC_COMM_SELF,MOFEM_DATA_INCONSISTENCT,"data inconsistency");
	    } 
	    if(nb_cols != transNN.size1()) {
	      SETERRQ(PETSC_COMM_SELF,MOFEM_DATA_INCONSISTENCT,"data inconsistency");
	    } 
	    ierr = MatSetValues(getFEMethod()->snes_B,nb_cols,cols,nb_rows,rows,trans_data,ADD_VALUES); CHKERRQ(ierr);
	  }

	}

	PetscFunctionReturn(0);
      }
    };

    /** \brief Rhs operaetar used to build matrix
      */
    template<typename FUNEVAL>
    struct OpRhs:public TriElementForcesAndSurcesCore::UserDataOperator {

      ublas::matrix<double> &hoCoords;
      FUNEVAL &functionEvaluator;
      int fieldNumber;

      OpRhs(const string field_name,ublas::matrix<double> &ho_coords,FUNEVAL &function_evaluator,int field_number): 
	TriElementForcesAndSurcesCore::UserDataOperator(field_name),
	hoCoords(ho_coords),functionEvaluator(function_evaluator),
	fieldNumber(field_number)  {}

      ublas::vector<FieldData> NTf;
      ublas::vector<DofIdx> iNdices;

      PetscErrorCode doWork(
	int side,EntityType type,DataForcesAndSurcesCore::EntData &data) {
	PetscFunctionBegin;
	PetscErrorCode ierr;
  
	try {

	  unsigned int nb_row = data.getIndices().size();
	  if(nb_row==0) PetscFunctionReturn(0);

	  const FENumeredDofMoFEMEntity *dof_ptr;
	  ierr = getMoFEMFEPtr()->get_row_dofs_by_petsc_gloabl_dof_idx(data.getIndices()[0],&dof_ptr); CHKERRQ(ierr);
	  unsigned int rank = dof_ptr->get_max_rank();
  
	  NTf.resize(nb_row/rank);
          iNdices.resize(nb_row/rank);

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
	    
	    ublas::vector<double>& a = functionEvaluator(x,y,z)[fieldNumber];

	    if(a.size()!=rank) {
	      SETERRQ(PETSC_COMM_SELF,MOFEM_DATA_INCONSISTENCT,"data inconsistency");
	    }


	    for(unsigned int rr = 0;rr<rank;rr++) {

	      ublas::noalias(iNdices) = ublas::vector_slice<ublas::vector<int> >
		(data.getIndices(), ublas::slice(rr, rank, data.getIndices().size()/rank));

	      noalias(NTf) = data.getN(gg,nb_row/rank)*a[rr]*val;
	      ierr = VecSetValues(getFEMethod()->snes_f,iNdices.size(),
		&iNdices[0],&*NTf.data().begin(),ADD_VALUES); CHKERRQ(ierr);

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

  struct DirichletBC : public DisplacementBCFEMethodPreAndPostProc {

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
	  SETERRQ(PETSC_COMM_SELF,1,"need to initialised from AnalyticalDirihletBC::solveProblem");
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

  template<typename FUNEVAL>
  PetscErrorCode setApproxOps(
    FieldInterface &m_field,
    string field_name,
    FUNEVAL &funtcion_evaluator,int field_number = 0,
    string nodals_positions = "MESH_NODE_POSITIONS") {
    PetscFunctionBegin;
    if(m_field.check_field(nodals_positions)) {
	approxField.getLoopFeApprox().get_op_to_do_Rhs().push_back(new ApproxField::OpHoCoord(nodals_positions,approxField.hoCoords));
    }
    approxField.getLoopFeApprox().get_op_to_do_Rhs().push_back(new ApproxField::OpRhs<FUNEVAL>(field_name,approxField.hoCoords,funtcion_evaluator,field_number));
    approxField.getLoopFeApprox().get_op_to_do_Lhs().push_back(new ApproxField::OpLhs(field_name,approxField.hoCoords));
    PetscFunctionReturn(0);
  }

  PetscErrorCode initializeProblem(
    FieldInterface &m_field,string fe,string field,
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
    ierr = m_field.add_ents_to_finite_element_by_TRIs(tRis,fe); CHKERRQ(ierr);
    PetscFunctionReturn(0);
  }

  Mat A;
  Vec D,F;
  KSP kspSolver;
  PetscErrorCode setProblem(
    FieldInterface &m_field,string problem) {
    PetscFunctionBegin;
    PetscErrorCode ierr;

    ierr = m_field.VecCreateGhost(problem,ROW,&F); CHKERRQ(ierr);
    ierr = m_field.VecCreateGhost(problem,COL,&D); CHKERRQ(ierr);
    ierr = m_field.MatCreateMPIAIJWithArrays(problem,&A); CHKERRQ(ierr);

    ierr = KSPCreate(PETSC_COMM_WORLD,&kspSolver); CHKERRQ(ierr);
    ierr = KSPSetOperators(kspSolver,A,A); CHKERRQ(ierr);
    ierr = KSPSetFromOptions(kspSolver); CHKERRQ(ierr);

    PC pc;
    ierr = KSPGetPC(kspSolver,&pc); CHKERRQ(ierr);
    ierr = PCSetType(pc,PCLU); CHKERRQ(ierr);
    ierr = PCFactorSetMatSolverPackage(pc,MATSOLVERMUMPS); CHKERRQ(ierr);
    ierr = PCFactorSetUpMatSolverPackage(pc);  CHKERRQ(ierr);

    PetscFunctionReturn(0);
  }
  
  PetscErrorCode solveProblem(
    FieldInterface &m_field,string problem,string fe,DirichletBC &bc) {
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

    ierr = KSPSolve(kspSolver,F,D); CHKERRQ(ierr);
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
    ierr = KSPDestroy(&kspSolver); CHKERRQ(ierr);
    ierr = MatDestroy(&A); CHKERRQ(ierr);
    ierr = VecDestroy(&F); CHKERRQ(ierr);
    ierr = VecDestroy(&D); CHKERRQ(ierr);
    PetscFunctionReturn(0);
  }



};

#endif //__ANALYTICALDIRIHLETBC_HPP__

