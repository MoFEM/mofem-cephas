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

#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/matrix_proxy.hpp>
#include <boost/numeric/ublas/vector.hpp>
#include <boost/numeric/ublas/vector_proxy.hpp>
#include <boost/numeric/ublas/io.hpp>
#include <boost/ptr_container/ptr_vector.hpp>

using namespace boost::numeric;

namespace MoFEM {

template<typename FUNEVAL>
struct FieldApproximation {

  FieldInterface &mField;
  const string problemName;
  VolumeH1H1ElementForcesAndSurcesCore fe;

  FieldApproximation(
    FieldInterface &m_field,
    const string &problem_name):
    mField(m_field),problemName(problemName),
    fe(m_field) {}

  Mat A;
  Vec D,F;

  PetscErrorCode createMatrixAndVector() {
    PetscFunctionBegin;
    PetscErrorCode ierr;
    ierr = mField.VecCreateGhost(problemName,Row,&F); CHKERRQ(ierr);
    ierr = mField.VecCreateGhost(problemName,Col,&D); CHKERRQ(ierr);
    ierr = mField.MatCreateMPIAIJWithArrays(problemName,&A); CHKERRQ(ierr);
    PetscFunctionReturn(0);
  }
  PetscErrorCode deleteMatrixAndVector() {
    PetscFunctionBegin;
    PetscErrorCode ierr;
    ierr = VecDestroy(&D); CHKERRQ(ierr);
    ierr = VecDestroy(&F); CHKERRQ(ierr);
    ierr = MatDestroy(&A); CHKERRQ(ierr);
    PetscFunctionReturn(0);
  }

  struct OpApprox: public VolumeH1H1ElementForcesAndSurcesCore::UserDataOperator {

    Mat &A;
    Vec &F;

    OpApprox(const string &field_name,Mat _A,Vec _F):
      VolumeH1H1ElementForcesAndSurcesCore::UserDataOperator(field_name),
      A(_A),F(_F),functionEvaluator(function_evaluator) {}

    ublas::matrix<FieldData> NN;
    ublas::vector<FieldData> Nf;

    PetscErrorCode doWork(
      int row_side,int col_side,
      EntityType row_type,EntityType col_type,
      DataForcesAndSurcesCore::EntData &row_data,
      DataForcesAndSurcesCore::EntData &col_data) {
      PetscFunctionBegin;

      PetscErrorCode ierr;

      const FENumeredDofMoFEMEntity *dof_ptr;
      ierr = ptrFE->fe_ptr->get_row_dofs_by_petsc_gloabl_dof_idx(row_data.getIndices()[0],&dof_ptr); CHKERRQ(ierr);
      int rank = dof_ptr->get_max_rank();

      int nb_row_dofs = row_data.getN().size2();
      int nb_col_dofs = col_data.getN().size2();

      NN.resize(nb_row_dofs,nb_col_dofs);
      bzero(&*NN.data().begin(),nb_row_dofs*nb_col_dofs*sizeof(FieldData));

      for(unsigned int gg = 0;gg<row_data.getN().size1();gg++) {
	double val = ptrFE->vOlume*ptrFE->gaussPts(gg,4);
	cblas_dger(CblasRowMajor,
	    nb_row_dofs,nb_col_dofs,
	    val,&row_data.getN()(gg,0),1,&col_data.getN()(gg,0),1,
	    &*NN.data().begin(),nb_col_dofs);

      }

      for(int rr = 0;rr < rank; rr++) {
      
	if((row_data.getIndices().size()%rank)!=0) {
	  SETERRQ(PETSC_COMM_SELF,1,"data inconsistency");
	}

	if((col_data.getIndices().size()%rank)!=0) {
	    SETERRQ(PETSC_COMM_SELF,1,"data inconsistency");
	}


	int nb_rows;
	int nb_cols;
	int *rows;
	int *cols;
	double *data = &*NN.data().begin();

	if(rank > 1) {

	  ublas::vector<DofIdx> row_indices( ublas::vector_slice<ublas::vector<DofIdx> >(row_data.getIndices(), ublas::slice(rr, rank, row_data.getIndices().size())));
	  ublas::vector<DofIdx> col_indices( ublas::vector_slice<ublas::vector<DofIdx> >(col_data.getIndices(), ublas::slice(rr, rank, col_data.getIndices().size())));

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

	ierr = MatSetValues(A,nb_rows,rows,nb_cols,cols,data,ADD_VALUES); CHKERRQ(ierr);

      }


      PetscFunctionReturn(0);	
    }



    PetscErrorCode doWork(
      int side,EntityType type,DataForcesAndSurcesCore::EntData &data) {
      PetscFunctionBegin;

      PetscErrorCode ierr;

      const FENumeredDofMoFEMEntity *dof_ptr;
      ierr = ptrFE->fe_ptr->get_row_dofs_by_petsc_gloabl_dof_idx(data.getIndices()[0],&dof_ptr); CHKERRQ(ierr);
      int rank = dof_ptr->get_max_rank();

      int nb_row_dofs = data.getN().size2();
      
      Nf.resize(data.getIndices().size());
    bzero(&*Nf.data().begin(),data.getIndices().size()*sizeof(double));

      for(unsigned int gg = 0;gg<data.getN().size1();gg++) {
      
	double x = ptrFE->coordsAtGaussPts(gg,0);
	double y = ptrFE->coordsAtGaussPts(gg,1);
	double z = ptrFE->coordsAtGaussPts(gg,2);
	ublas::vector<FieldData> fun_val = functionEvaluator(x,y,z);
	double val = ptrFE->vOlume*ptrFE->gaussPts(gg,4);
	
	if(fun_val.size() != rank) {
	  SETERRQ(PETSC_COMM_SELF,1,"data inconsistency");
	}
	for(int rr = 0;rr<rank;rr++) {
	  cblas_daxpy(nb_row_dofs,val*fun_val[rr],&data.getN()(gg,0),1,&Nf[rr],rank);
	}

      }

      for(int rr = 0;rr < rank; rr++) {

	int nb_rows;
	int *rows;
	double *f_data = &*Nf.data().begin();

	if(rank > 1) {

	  ublas::vector<DofIdx> indices( ublas::vector_slice<ublas::vector<DofIdx> >(data.getIndices(), ublas::slice(rr, rank, data.getIndices().size())));

	  nb_rows = indices.size();
	  rows = &*indices.data().begin();

	} else {

	  nb_rows = data.getIndices().size();
	  rows = &*data.getIndices().data().begin();

	}

	ierr = VecSetValues(F,nb_rows,rows,f_data,ADD_VALUES); CHKERRQ(ierr);
      
      }

      PetscFunctionReturn(0);
    }


  };

  PetscErrorCode computeMatrixAndVector() {
    PetscFunctionBegin;

    OpApprox approx(my_split);
    fe1.get_op_to_do_NH1().push_back(&op);
    fe1.get_op_to_do_NH1NH1().push_back(&op);


    PetscFunctionReturn(0);
  }
  PetscErrorCode solevProblem() {
    PetscFunctionBegin;
    PetscFunctionReturn(0);
  }

};


}

#endif //__FILEDAPPROXIMATION_HPP

