/* Copyright (C) 2013, Lukasz Kaczmarczyk (likask AT wp.pl)
* --------------------------------------------------------------
*
* DESCRIPTION: FIXME
*
* This is not exactly procedure for linear elatic dynamics, since jacobian is
* evaluated at every time step and snes procedure is involved. However it is
* implemented like that, to test methodology for general nonlinear problem.
*
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

#ifndef __FILEDAPPROXIMATION_HPP
#define __FILEDAPPROXIMATION_HPP

using namespace boost::numeric;

/** \brief Finite element for approximating analytical filed on the mesh
  * \ingroup user_modules
  */
template<typename FUNEVAL>
struct FieldApproximationH1 {

  FieldInterface &mField;
  const string problemName;
  TetElementForcesAndSourcesCore fe;

  FieldApproximationH1(
    FieldInterface &m_field):
    mField(m_field),fe(m_field) {}

  /** \brief Gauss point poperatiors to caculete matrices and vectors
    *
    * Function work on thetrahedrals
    */
  struct OpApprox: public TetElementForcesAndSourcesCore::UserDataOperator {

    Mat A;
    Vec F;
    FUNEVAL &functionEvaluator;

    OpApprox(const string &field_name,Mat _A,Vec _F,FUNEVAL &function_evaluator):
      TetElementForcesAndSourcesCore::UserDataOperator(field_name),
      A(_A),F(_F),functionEvaluator(function_evaluator) {}
    ~OpApprox() {}

    ublas::matrix<FieldData> NN;
    ublas::matrix<FieldData> transNN;
    ublas::vector<FieldData> Nf;

    /** \brief calulate matrix
      */
    PetscErrorCode doWork(
      int row_side,int col_side,
      EntityType row_type,EntityType col_type,
      DataForcesAndSurcesCore::EntData &row_data,
      DataForcesAndSurcesCore::EntData &col_data) {
      PetscFunctionBegin;


	  //std::string wait;
	  //std::cout << "\n I am here !!!!!!! rank = \n" << std::endl;

      if(row_data.getIndices().size()==0) PetscFunctionReturn(0);
      if(col_data.getIndices().size()==0) PetscFunctionReturn(0);

      PetscErrorCode ierr;

      const FENumeredDofMoFEMEntity *dof_ptr;
      ierr = getMoFEMFEPtr()->get_row_dofs_by_petsc_gloabl_dof_idx(row_data.getIndices()[0],&dof_ptr); CHKERRQ(ierr);
      int rank = dof_ptr->get_max_rank();

      int nb_row_dofs = row_data.getIndices().size()/rank;
      int nb_col_dofs = col_data.getIndices().size()/rank;

      NN.resize(nb_row_dofs,nb_col_dofs);
      bzero(&*NN.data().begin(),nb_row_dofs*nb_col_dofs*sizeof(FieldData));





      unsigned int nb_gauss_pts = row_data.getN().size1();
      for(unsigned int gg = 0;gg<nb_gauss_pts;gg++) {
	double w = getVolume()*getGaussPts()(3,gg);
	if(getHoCoordsAtGaussPts().size1()==nb_gauss_pts) {
	  w *= getHoGaussPtsDetJac()[gg];
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
	  SETERRQ(PETSC_COMM_SELF,1,"data inconsistency");
	}

	if((col_data.getIndices().size()%rank)!=0) {
	    SETERRQ(PETSC_COMM_SELF,1,"data inconsistency");
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
	  SETERRQ(PETSC_COMM_SELF,1,"data inconsistency");
	}
	if(nb_cols != NN.size2()) {
	  SETERRQ(PETSC_COMM_SELF,1,"data inconsistency");
	}

	ierr = MatSetValues(A,nb_rows,rows,nb_cols,cols,data,ADD_VALUES); CHKERRQ(ierr);
	if( (row_type != col_type) || (row_side != col_side) ) {
	  if(nb_rows != transNN.size2()) {
	    SETERRQ(PETSC_COMM_SELF,1,"data inconsistency");
	  }
	  if(nb_cols != transNN.size1()) {
	    SETERRQ(PETSC_COMM_SELF,1,"data inconsistency");
	  }
	  ierr = MatSetValues(A,nb_cols,cols,nb_rows,rows,trans_data,ADD_VALUES); CHKERRQ(ierr);
	}

      }

      PetscFunctionReturn(0);
    }

    /** \brief caclulate vector
      */
    PetscErrorCode doWork(
      int side,EntityType type,DataForcesAndSurcesCore::EntData &data) {
      PetscFunctionBegin;

      if(data.getIndices().size()==0) PetscFunctionReturn(0);

      PetscErrorCode ierr;

      //PetscAttachDebugger();

      const FENumeredDofMoFEMEntity *dof_ptr;
      ierr = getMoFEMFEPtr()->get_row_dofs_by_petsc_gloabl_dof_idx(data.getIndices()[0],&dof_ptr); CHKERRQ(ierr);
      unsigned int rank = dof_ptr->get_max_rank();

      int nb_row_dofs = data.getIndices().size()/rank;

      Nf.resize(data.getIndices().size());
      bzero(&*Nf.data().begin(),data.getIndices().size()*sizeof(FieldData));

      if(getCoordsAtGaussPts().size1()!=data.getN().size1()) {
	SETERRQ(PETSC_COMM_SELF,1,"data inconsistency");
      }
      if(getCoordsAtGaussPts().size2()!=3) {
	SETERRQ(PETSC_COMM_SELF,1,"data inconsistency");
      }

      // itegration
      unsigned int nb_gauss_pts = data.getN().size1();
      for(unsigned int gg = 0;gg<nb_gauss_pts;gg++) {

	double x,y,z,w;
	w = getVolume()*getGaussPts()(3,gg);
	if(getHoCoordsAtGaussPts().size1()==nb_gauss_pts) {
	  //intergation poits global positions if higher order geometry is given
	  x = getHoCoordsAtGaussPts()(gg,0);
	  y = getHoCoordsAtGaussPts()(gg,1);
	  z = getHoCoordsAtGaussPts()(gg,2);
	  // correction of jacobian for higher order geometry
	  w *= getHoGaussPtsDetJac()[gg];
	} else {
	  //intergartion point global positions for linear tetrahedral element
	  x = getCoordsAtGaussPts()(gg,0);
	  y = getCoordsAtGaussPts()(gg,1);
	  z = getCoordsAtGaussPts()(gg,2);
	}

	//cerr << x << " " << y << " " << z << " " << w << " " << getHoGaussPtsDetJac()[gg] << endl;

	ublas::vector<FieldData> fun_val = functionEvaluator(x,y,z);

	if(fun_val.size() != rank) {
	  SETERRQ(PETSC_COMM_SELF,1,"data inconsistency");
	}

	for(unsigned int rr = 0;rr<rank;rr++) {
	  cblas_daxpy(nb_row_dofs,w*fun_val[rr],&data.getN()(gg,0),1,&Nf[rr],rank);
	}

      }

      ierr = VecSetValues(F,data.getIndices().size(),
	&data.getIndices()[0],&Nf[0],ADD_VALUES); CHKERRQ(ierr);

      PetscFunctionReturn(0);
    }


  };

  /** \brief assemble matrix and vector
    */
  PetscErrorCode loopMatrixAndVector(
    const string &problem_name,const string &fe_name,const string &field_name,Mat A,Vec F,
    FUNEVAL &function_evaluator) {
    PetscFunctionBegin;
    PetscErrorCode ierr;

    //add operator to calulate F vector
    fe.get_op_to_do_Rhs().push_back(new OpApprox(field_name,A,F,function_evaluator));
    //add operator to calulate A matrix
    fe.get_op_to_do_Lhs().push_back(new OpApprox(field_name,A,F,function_evaluator));

    MatZeroEntries(A);
    VecZeroEntries(F);
    ierr = VecGhostUpdateBegin(F,INSERT_VALUES,SCATTER_FORWARD); CHKERRQ(ierr);
    ierr = VecGhostUpdateEnd(F,INSERT_VALUES,SCATTER_FORWARD); CHKERRQ(ierr);
    //calulate and assembe
    ierr = mField.loop_finite_elements(problem_name,fe_name,fe);  CHKERRQ(ierr);
    ierr = MatAssemblyBegin(A,MAT_FLUSH_ASSEMBLY); CHKERRQ(ierr);
    ierr = MatAssemblyEnd(A,MAT_FLUSH_ASSEMBLY); CHKERRQ(ierr);
    ierr = VecAssemblyBegin(F); CHKERRQ(ierr);
    ierr = VecAssemblyEnd(F); CHKERRQ(ierr);
    PetscFunctionReturn(0);
  }

};

#endif //__FILEDAPPROXIMATION_HPP
