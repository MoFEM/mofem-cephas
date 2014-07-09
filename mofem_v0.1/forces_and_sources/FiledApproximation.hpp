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

#include "FieldInterface.hpp"
#include "FEM.h"

#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/matrix_proxy.hpp>
#include <boost/numeric/ublas/vector.hpp>
#include <boost/numeric/ublas/vector_proxy.hpp>
#include <boost/numeric/ublas/io.hpp>
#include <boost/ptr_container/ptr_vector.hpp>

using namespace boost::numeric;

namespace MoFEM {

template<typename FUNEVAL>
struct FieldApproximationH1 {

  FieldInterface &mField;
  const string problemName;
  TetElementForcesAndSourcesCore fe;

  FieldApproximationH1(
    FieldInterface &m_field):
    mField(m_field),fe(m_field) {}

  struct OpApprox: public TetElementForcesAndSourcesCore::UserDataOperator {

    Mat &A;
    Vec &F;
    FUNEVAL &functionEvaluator;

    OpApprox(const string &field_name,Mat &_A,Vec &_F,FUNEVAL &function_evaluator):
      TetElementForcesAndSourcesCore::UserDataOperator(field_name),
      A(_A),F(_F),functionEvaluator(function_evaluator) {}
    ~OpApprox() {}

    ublas::matrix<FieldData> NN;
    ublas::matrix<FieldData> transNN;
    ublas::vector<FieldData> Nf;

    PetscErrorCode doWork(
      int row_side,int col_side,
      EntityType row_type,EntityType col_type,
      DataForcesAndSurcesCore::EntData &row_data,
      DataForcesAndSurcesCore::EntData &col_data) {
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
      bzero(&*NN.data().begin(),nb_row_dofs*nb_col_dofs*sizeof(FieldData));
      
      for(unsigned int gg = 0;gg<row_data.getN().size1();gg++) {
	double val = getVolume()*getGaussPts()(3,gg);
	cblas_dger(CblasRowMajor,
	    nb_row_dofs,nb_col_dofs,
	    val,&row_data.getN()(gg,0),1,&col_data.getN()(gg,0),1,
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

      for(unsigned int gg = 0;gg<data.getN().size1();gg++) {
      
	double x = getCoordsAtGaussPts()(gg,0);
	double y = getCoordsAtGaussPts()(gg,1);
	double z = getCoordsAtGaussPts()(gg,2);
	ublas::vector<FieldData> fun_val = functionEvaluator(x,y,z);
	if(fun_val.size() != rank) {
	  SETERRQ(PETSC_COMM_SELF,1,"data inconsistency");
	}

	double val = getVolume()*getGaussPts()(3,gg);
	for(unsigned int rr = 0;rr<rank;rr++) {
	  cblas_daxpy(nb_row_dofs,val*fun_val[rr],&data.getN()(gg,0),1,&Nf[rr],rank);
	}

      }

      ierr = VecSetValues(F,data.getIndices().size(),
	&data.getIndices()[0],&Nf[0],ADD_VALUES); CHKERRQ(ierr);

      PetscFunctionReturn(0);
    }


  };

  PetscErrorCode loopMatrixAndVector(
    const string &problem_name,const string &fe_name,const string &field_name,Mat A,Vec F,
    FUNEVAL &function_evaluator) {
    PetscFunctionBegin;
    PetscErrorCode ierr;
    fe.get_op_to_do_Rhs_H1().push_back(new OpApprox(field_name,A,F,function_evaluator));
    fe.get_op_to_do_Lhs_H1H1().push_back(new OpApprox(field_name,A,F,function_evaluator));
    MatZeroEntries(A);
    VecZeroEntries(F);
    ierr = VecGhostUpdateBegin(F,INSERT_VALUES,SCATTER_FORWARD); CHKERRQ(ierr);
    ierr = VecGhostUpdateEnd(F,INSERT_VALUES,SCATTER_FORWARD); CHKERRQ(ierr);
    ierr = mField.loop_finiteElementsPtr(problem_name,"TEST_FE",fe);  CHKERRQ(ierr);
    ierr = MatAssemblyBegin(A,MAT_FLUSH_ASSEMBLY); CHKERRQ(ierr);
    ierr = MatAssemblyEnd(A,MAT_FLUSH_ASSEMBLY); CHKERRQ(ierr);
    ierr = VecAssemblyBegin(F); CHKERRQ(ierr);
    ierr = VecAssemblyEnd(F); CHKERRQ(ierr);
    PetscFunctionReturn(0);
  }

};


}

#endif //__FILEDAPPROXIMATION_HPP

