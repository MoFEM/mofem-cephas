/* \file FieldApproximation.hpp 

\brief Element to calculate approximation on volume elements

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

namespace MoFEM {

/** \brief Finite element for approximating analytical filed on the mesh
  * \ingroup mofem_forces_and_sources
  */
struct FieldApproximationH1 {

  FieldInterface &mField;
  const string problemName;
  VolumeElementForcesAndSourcesCore fe;
  int addToRule;

  FieldApproximationH1(
    FieldInterface &m_field):
    mField(m_field),fe(m_field),addToRule(0) {}

  /** \brief set integration rule

    Note: Add 1 to take into account 2nd order geometry approximation. Overload
    that function if linear or HO approximation is set.

  */
  int getRule(int order) { return order+addToRule; }; 

  /** \brief Gauss point poperatiors to caculete matrices and vectors
    *
    * Function work on thetrahedrals
    */
  template<typename FUNEVAL>
  struct OpApprox: public VolumeElementForcesAndSourcesCore::UserDataOperator {

    Mat A;
    vector<Vec> &vecF;
    FUNEVAL &functionEvaluator;

    OpApprox(const string &field_name,Mat _A,vector<Vec> &vec_F,FUNEVAL &function_evaluator):
      VolumeElementForcesAndSourcesCore::UserDataOperator(field_name),
      A(_A),vecF(vec_F),functionEvaluator(function_evaluator) {}
    virtual ~OpApprox() {}

    ublas::matrix<FieldData> NN;
    ublas::matrix<FieldData> transNN;
    vector<ublas::vector<FieldData> > Nf;

    /** \brief calulate matrix
      */
    PetscErrorCode doWork(
      int row_side,int col_side,
      EntityType row_type,EntityType col_type,
      DataForcesAndSurcesCore::EntData &row_data,
      DataForcesAndSurcesCore::EntData &col_data) {
      PetscFunctionBegin;

      if(A == PETSC_NULL) PetscFunctionReturn(0);
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

	ierr = MatSetValues(A,nb_rows,rows,nb_cols,cols,data,ADD_VALUES); CHKERRQ(ierr);
	if( (row_type != col_type) || (row_side != col_side) ) {
	  if(nb_rows != transNN.size2()) {
	    SETERRQ(PETSC_COMM_SELF,MOFEM_DATA_INCONSISTENCT,"data inconsistency");
	  } 
	  if(nb_cols != transNN.size1()) {
	    SETERRQ(PETSC_COMM_SELF,MOFEM_DATA_INCONSISTENCT,"data inconsistency");
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
      
      if(getCoordsAtGaussPts().size1()!=data.getN().size1()) {
	SETERRQ(PETSC_COMM_SELF,MOFEM_DATA_INCONSISTENCT,"data inconsistency");
      }
      if(getCoordsAtGaussPts().size2()!=3) {
	SETERRQ(PETSC_COMM_SELF,MOFEM_DATA_INCONSISTENCT,"data inconsistency");
      }

      // itegration 
      unsigned int nb_gauss_pts = data.getN().size1();
      for(unsigned int gg = 0;gg != nb_gauss_pts;gg++) {

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

	vector<ublas::vector<FieldData> > fun_val;
	try {

	   fun_val = functionEvaluator(x,y,z);

	} catch (exception& ex) {
	  ostringstream ss;
	  ss << "thorw in method: " << ex.what() << " at line " << __LINE__ << " in file " << __FILE__;
	  SETERRQ(PETSC_COMM_SELF,MOFEM_STD_EXCEPTION_THROW,ss.str().c_str());
	}

	if(fun_val.size()!=vecF.size()) {
	  SETERRQ(PETSC_COMM_SELF,MOFEM_DATA_INCONSISTENCT,"data inconsistency");

	}

	Nf.resize(fun_val.size());
	for(unsigned int lhs = 0;lhs != fun_val.size();lhs++) {
	  
	  if(!gg) {
	    Nf[lhs].resize(data.getIndices().size());
	    Nf[lhs].clear();
	  } 
		
	  if(fun_val[lhs].size() != rank) {
	    SETERRQ(PETSC_COMM_SELF,MOFEM_DATA_INCONSISTENCT,"data inconsistency");
	  }

	  for(unsigned int rr = 0;rr != rank;rr++) {
	    cblas_daxpy(nb_row_dofs,w*(fun_val[lhs])[rr],&data.getN()(gg,0),1,&(Nf[lhs])[rr],rank);
	  }

	}

      }

      for(unsigned int lhs = 0;lhs != vecF.size();lhs++) {

	ierr = VecSetValues(vecF[lhs],data.getIndices().size(),
	  &data.getIndices()[0],&(Nf[lhs])[0],ADD_VALUES); CHKERRQ(ierr);

      }

      PetscFunctionReturn(0);
    }


  };

  /** \brief assemble matrix and vector 
    */
  template<typename FUNEVAL>
  PetscErrorCode loopMatrixAndVector(
    const string &problem_name,const string &fe_name,const string &field_name,
    Mat A,vector<Vec> &vec_F,FUNEVAL &function_evaluator) {
    PetscFunctionBegin;
    PetscErrorCode ierr;
	
    //add operator to calulate F vector
    fe.get_op_to_do_Rhs().push_back(new OpApprox<FUNEVAL>(field_name,A,vec_F,function_evaluator));
    //add operator to calulate A matrix
    if(A) {
      fe.get_op_to_do_Lhs().push_back(new OpApprox<FUNEVAL>(field_name,A,vec_F,function_evaluator));
    }
	
    if(A) {
      ierr = MatZeroEntries(A); CHKERRQ(ierr);
    }

    for(unsigned int lhs = 0; lhs<vec_F.size(); lhs++) {

      ierr = VecZeroEntries(vec_F[lhs]); CHKERRQ(ierr);
      ierr = VecGhostUpdateBegin(vec_F[lhs],INSERT_VALUES,SCATTER_FORWARD); CHKERRQ(ierr);
      ierr = VecGhostUpdateEnd(vec_F[lhs],INSERT_VALUES,SCATTER_FORWARD); CHKERRQ(ierr);

    }
    
    //calulate and assembe
    ierr = mField.loop_finite_elements(problem_name,fe_name,fe);  CHKERRQ(ierr); 
    
    if(A) {
      ierr = MatAssemblyBegin(A,MAT_FLUSH_ASSEMBLY); CHKERRQ(ierr);
      ierr = MatAssemblyEnd(A,MAT_FLUSH_ASSEMBLY); CHKERRQ(ierr);
    }

    for(unsigned int lhs = 0; lhs<vec_F.size(); lhs++) {

      ierr = VecAssemblyBegin(vec_F[lhs]); CHKERRQ(ierr);
      ierr = VecAssemblyEnd(vec_F[lhs]); CHKERRQ(ierr);

    }
  
    PetscFunctionReturn(0);
  }

};


}

#endif //__FILEDAPPROXIMATION_HPP

