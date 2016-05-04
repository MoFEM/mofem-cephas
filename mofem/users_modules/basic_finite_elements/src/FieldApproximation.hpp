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

/** \brief Finite element for approximating analytical filed on the mesh
  * \ingroup user_modules
  */
struct FieldApproximationH1 {

  FieldInterface &mField;
  const std::string problemName;
  VolumeElementForcesAndSourcesCore feVolume;
  FaceElementForcesAndSourcesCore feFace;

  FieldApproximationH1(
    FieldInterface &m_field):
    mField(m_field),
    feVolume(m_field),
    feFace(m_field) {
  }

  /** \brief Gauss point operators to calculate matrices and vectors
    *
    * Function work on volumes (Terahedrons & Bricks)
    */
  template<typename FUNEVAL>
  struct OpApproxVolume: public VolumeElementForcesAndSourcesCore::UserDataOperator {

    Mat A;
    std::vector<Vec> &vecF;
    FUNEVAL &functionEvaluator;

    OpApproxVolume(const std::string &field_name,Mat _A,std::vector<Vec> &vec_F,FUNEVAL &function_evaluator):
      VolumeElementForcesAndSourcesCore::UserDataOperator(field_name,UserDataOperator::OPROW|UserDataOperator::OPROWCOL),
      A(_A),vecF(vec_F),functionEvaluator(function_evaluator) {}
    virtual ~OpApproxVolume() {}

    MatrixDouble NN;
    MatrixDouble transNN;
    std::vector<VectorDouble > Nf;

    /** \brief calculate matrix
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

      const FENumeredDofEntity *dof_ptr;
      ierr = getMoFEMFEPtr()->get_row_dofs_by_petsc_gloabl_dof_idx(row_data.getIndices()[0],&dof_ptr); CHKERRQ(ierr);
      int rank = dof_ptr->get_nb_of_coeffs();

      int nb_row_dofs = row_data.getIndices().size()/rank;
      int nb_col_dofs = col_data.getIndices().size()/rank;

      NN.resize(nb_row_dofs,nb_col_dofs,false);
      NN.clear();

      unsigned int nb_gauss_pts = row_data.getN().size1();
      for(unsigned int gg = 0;gg!=nb_gauss_pts;gg++) {
        double w = getVolume()*getGaussPts()(3,gg);
        if(getHoCoordsAtGaussPts().size1()==nb_gauss_pts) {
          w *= getHoGaussPtsDetJac()[gg];
        }
        //noalias(NN) += w*outer_prod(row_data.getN(gg),col_data.getN(gg));
        cblas_dger(CblasRowMajor,
          nb_row_dofs,nb_col_dofs,
          w,&row_data.getN()(gg,0),1,&col_data.getN()(gg,0),1,
          &*NN.data().begin(),nb_col_dofs
        );
      }

      if( (row_type != col_type) || (row_side != col_side) ) {
        transNN.resize(nb_col_dofs,nb_row_dofs,false);
        ublas::noalias(transNN) = trans(NN);
      }

      double *data = &*NN.data().begin();
      double *trans_data = &*transNN.data().begin();
      VectorInt row_indices,col_indices;
      row_indices.resize(nb_row_dofs);
      col_indices.resize(nb_col_dofs);

      for(int rr = 0;rr < rank; rr++) {

        if((row_data.getIndices().size()%rank)!=0) {
          SETERRQ(PETSC_COMM_SELF,MOFEM_DATA_INCONSISTENCY,"data inconsistency");
        }

        if((col_data.getIndices().size()%rank)!=0) {
          SETERRQ(PETSC_COMM_SELF,MOFEM_DATA_INCONSISTENCY,"data inconsistency");
        }

        unsigned int nb_rows;
        unsigned int nb_cols;
        int *rows;
        int *cols;

        if(rank > 1) {

          ublas::noalias(row_indices) = ublas::vector_slice<VectorInt >
          (row_data.getIndices(), ublas::slice(rr, rank, row_data.getIndices().size()/rank));
          ublas::noalias(col_indices) = ublas::vector_slice<VectorInt >
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
          SETERRQ(PETSC_COMM_SELF,MOFEM_DATA_INCONSISTENCY,"data inconsistency");
        }
        if(nb_cols != NN.size2()) {
          SETERRQ(PETSC_COMM_SELF,MOFEM_DATA_INCONSISTENCY,"data inconsistency");
        }

        ierr = MatSetValues(A,nb_rows,rows,nb_cols,cols,data,ADD_VALUES); CHKERRQ(ierr);
        if( (row_type != col_type) || (row_side != col_side) ) {
          if(nb_rows != transNN.size2()) {
            SETERRQ(PETSC_COMM_SELF,MOFEM_DATA_INCONSISTENCY,"data inconsistency");
          }
          if(nb_cols != transNN.size1()) {
            SETERRQ(PETSC_COMM_SELF,MOFEM_DATA_INCONSISTENCY,"data inconsistency");
          }
          ierr = MatSetValues(A,nb_cols,cols,nb_rows,rows,trans_data,ADD_VALUES); CHKERRQ(ierr);
        }

      }

      PetscFunctionReturn(0);
    }

    /** \brief calculate vector
      */
    PetscErrorCode doWork(
      int side,EntityType type,DataForcesAndSurcesCore::EntData &data) {
      PetscFunctionBegin;

      if(data.getIndices().size()==0) PetscFunctionReturn(0);

      PetscErrorCode ierr;

      //PetscAttachDebugger();

      const FENumeredDofEntity *dof_ptr;
      ierr = getMoFEMFEPtr()->get_row_dofs_by_petsc_gloabl_dof_idx(data.getIndices()[0],&dof_ptr); CHKERRQ(ierr);
      unsigned int rank = dof_ptr->get_nb_of_coeffs();

      int nb_row_dofs = data.getIndices().size()/rank;

      if(getCoordsAtGaussPts().size1()!=data.getN().size1()) {
        SETERRQ(PETSC_COMM_SELF,MOFEM_DATA_INCONSISTENCY,"data inconsistency");
      }
      if(getCoordsAtGaussPts().size2()!=3) {
        SETERRQ(PETSC_COMM_SELF,MOFEM_DATA_INCONSISTENCY,"data inconsistency");
      }

      // itegration
      unsigned int nb_gauss_pts = data.getN().size1();
      for(unsigned int gg = 0;gg != nb_gauss_pts;gg++) {

        double x,y,z,w;
        w = getVolume()*getGaussPts()(3,gg);
        if(getHoCoordsAtGaussPts().size1()==nb_gauss_pts) {
          //intergation points global positions if higher order geometry is given
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

        std::vector<ublas::vector<double> > fun_val;
        try {

          fun_val = functionEvaluator(x,y,z);

        } catch (std::exception& ex) {
          std::ostringstream ss;
          ss << "thorw in method: " << ex.what() << " at line " << __LINE__ << " in file " << __FILE__;
          SETERRQ(PETSC_COMM_SELF,MOFEM_STD_EXCEPTION_THROW,ss.str().c_str());
        }

        if(fun_val.size()!=vecF.size()) {
          SETERRQ(PETSC_COMM_SELF,MOFEM_DATA_INCONSISTENCY,"data inconsistency");

        }

        Nf.resize(fun_val.size());
        for(unsigned int lhs = 0;lhs != fun_val.size();lhs++) {

          if(!gg) {
            Nf[lhs].resize(data.getIndices().size());
            Nf[lhs].clear();
          }

          if(fun_val[lhs].size() != rank) {
            SETERRQ(PETSC_COMM_SELF,MOFEM_DATA_INCONSISTENCY,"data inconsistency");
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

  /** \brief Gauss point operators to calculate matrices and vectors
    *
    * Function work on faces (Triangles & Quads)
    */
  template<typename FUNEVAL>
  struct OpApproxFace: public FaceElementForcesAndSourcesCore::UserDataOperator {

    Mat A;
    std::vector<Vec> &vecF;
    FUNEVAL &functionEvaluator;

    OpApproxFace(const std::string &field_name,Mat _A,std::vector<Vec> &vec_F,FUNEVAL &function_evaluator):
      FaceElementForcesAndSourcesCore::UserDataOperator(
        field_name,UserDataOperator::OPROW|UserDataOperator::OPROWCOL
      ),
      A(_A),
      vecF(vec_F),
      functionEvaluator(function_evaluator) {}
    virtual ~OpApproxFace() {}

    ublas::matrix<double> NN;
    ublas::matrix<double> transNN;
    std::vector<ublas::vector<double> > Nf;

    /** \brief calculate matrix
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
      const FENumeredDofEntity *dof_ptr;
      ierr = getMoFEMFEPtr()->get_row_dofs_by_petsc_gloabl_dof_idx(row_data.getIndices()[0],&dof_ptr); CHKERRQ(ierr);
      int rank = dof_ptr->get_nb_of_coeffs();
      int nb_row_dofs = row_data.getIndices().size()/rank;
      int nb_col_dofs = col_data.getIndices().size()/rank;
      NN.resize(nb_row_dofs,nb_col_dofs,false);
      NN.clear();
      unsigned int nb_gauss_pts = row_data.getN().size1();
      for(unsigned int gg = 0;gg!=nb_gauss_pts;gg++) {
        double w = getGaussPts()(2,gg);
        if(getNormals_at_GaussPt().size1()) {
          w *= 0.5*cblas_dnrm2(3,&getNormals_at_GaussPt()(gg,0),1);
        } else {
          w *= getArea();
        }
        //noalias(NN) += w*outer_prod(row_data.getN(gg),col_data.getN(gg));
        cblas_dger(CblasRowMajor,
          nb_row_dofs,nb_col_dofs,
          w,&row_data.getN()(gg,0),1,&col_data.getN()(gg,0),1,
          &*NN.data().begin(),nb_col_dofs
        );
      }

      if( (row_type != col_type) || (row_side != col_side) ) {
        transNN.resize(nb_col_dofs,nb_row_dofs,false);
        ublas::noalias(transNN) = trans(NN);
      }

      double *data = &*NN.data().begin();
      double *trans_data = &*transNN.data().begin();
      VectorInt row_indices,col_indices;
      row_indices.resize(nb_row_dofs);
      col_indices.resize(nb_col_dofs);
      for(int rr = 0;rr < rank; rr++) {
        if((row_data.getIndices().size()%rank)!=0) {
          SETERRQ(PETSC_COMM_SELF,MOFEM_DATA_INCONSISTENCY,"data inconsistency");
        }
        if((col_data.getIndices().size()%rank)!=0) {
          SETERRQ(PETSC_COMM_SELF,MOFEM_DATA_INCONSISTENCY,"data inconsistency");
        }
        unsigned int nb_rows;
        unsigned int nb_cols;
        int *rows;
        int *cols;
        if(rank > 1) {
          ublas::noalias(row_indices) = ublas::vector_slice<VectorInt >
          (row_data.getIndices(), ublas::slice(rr, rank, row_data.getIndices().size()/rank));
          ublas::noalias(col_indices) = ublas::vector_slice<VectorInt >
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
          SETERRQ(PETSC_COMM_SELF,MOFEM_DATA_INCONSISTENCY,"data inconsistency");
        }
        if(nb_cols != NN.size2()) {
          SETERRQ(PETSC_COMM_SELF,MOFEM_DATA_INCONSISTENCY,"data inconsistency");
        }
        ierr = MatSetValues(A,nb_rows,rows,nb_cols,cols,data,ADD_VALUES); CHKERRQ(ierr);
        if( (row_type != col_type) || (row_side != col_side) ) {
          if(nb_rows != transNN.size2()) {
            SETERRQ(PETSC_COMM_SELF,MOFEM_DATA_INCONSISTENCY,"data inconsistency");
          }
          if(nb_cols != transNN.size1()) {
            SETERRQ(PETSC_COMM_SELF,MOFEM_DATA_INCONSISTENCY,"data inconsistency");
          }
          ierr = MatSetValues(A,nb_cols,cols,nb_rows,rows,trans_data,ADD_VALUES); CHKERRQ(ierr);
        }
      }
      PetscFunctionReturn(0);
    }

    /** \brief calculate vector
      */
    PetscErrorCode doWork(
      int side,EntityType type,DataForcesAndSurcesCore::EntData &data) {
      PetscFunctionBegin;

      if(data.getIndices().size()==0) PetscFunctionReturn(0);
      PetscErrorCode ierr;

      //PetscAttachDebugger();
      const FENumeredDofEntity *dof_ptr;
      ierr = getMoFEMFEPtr()->get_row_dofs_by_petsc_gloabl_dof_idx(data.getIndices()[0],&dof_ptr); CHKERRQ(ierr);
      unsigned int rank = dof_ptr->get_nb_of_coeffs();

      int nb_row_dofs = data.getIndices().size()/rank;

      if(getCoordsAtGaussPts().size1()!=data.getN().size1()) {
        SETERRQ(PETSC_COMM_SELF,MOFEM_DATA_INCONSISTENCY,"data inconsistency");
      }
      if(getCoordsAtGaussPts().size2()!=3) {
        SETERRQ(PETSC_COMM_SELF,MOFEM_DATA_INCONSISTENCY,"data inconsistency");
      }

      VectorDouble normal(3);
      VectorDouble tangent1(3);
      VectorDouble tangent2(3);
      tangent1.clear();
      tangent2.clear();

      // integration
      unsigned int nb_gauss_pts = data.getN().size1();
      for(unsigned int gg = 0;gg != nb_gauss_pts;gg++) {
        double x,y,z,w;
        w = getGaussPts()(2,gg);
        if(getNormals_at_GaussPt().size1()) {
          w *= 0.5*cblas_dnrm2(3,&getNormals_at_GaussPt()(gg,0),1);
        } else {
          w *= getArea();
        }

        if(getHoCoordsAtGaussPts().size1()==nb_gauss_pts) {
          //intergation points global positions if higher order geometry is given
          x = getHoCoordsAtGaussPts()(gg,0);
          y = getHoCoordsAtGaussPts()(gg,1);
          z = getHoCoordsAtGaussPts()(gg,2);
        } else {
          //intergartion point global positions for linear tetrahedral element
          x = getCoordsAtGaussPts()(gg,0);
          y = getCoordsAtGaussPts()(gg,1);
          z = getCoordsAtGaussPts()(gg,2);
        }

        if(getNormals_at_GaussPt().size1()) {
         noalias(normal) = getNormals_at_GaussPt(gg);
         for(int dd = 0;dd<3;dd++) {
          tangent1[dd] = getTangent1_at_GaussPt()(gg,dd);
          tangent2[dd] = getTangent2_at_GaussPt()(gg,dd);
        }
       } else {
         noalias(normal) = getNormal();
       }

        std::vector<ublas::vector<double> > fun_val;
        try {
          EntityHandle ent = getFEMethod()->fePtr->get_ent();
          fun_val = functionEvaluator(ent,x,y,z,normal,tangent1,tangent2);
        } catch (std::exception& ex) {
          std::ostringstream ss;
          ss << "thorw in method: " << ex.what() << " at line " << __LINE__ << " in file " << __FILE__;
          SETERRQ(PETSC_COMM_SELF,MOFEM_STD_EXCEPTION_THROW,ss.str().c_str());
        }

        if(fun_val.size()!=vecF.size()) {
          SETERRQ(PETSC_COMM_SELF,MOFEM_DATA_INCONSISTENCY,"data inconsistency");
        }

        Nf.resize(fun_val.size());
        for(unsigned int lhs = 0;lhs != fun_val.size();lhs++) {
          if(!gg) {
            Nf[lhs].resize(data.getIndices().size());
            Nf[lhs].clear();
          }
          if(fun_val[lhs].size() != rank) {
            SETERRQ(PETSC_COMM_SELF,MOFEM_DATA_INCONSISTENCY,"data inconsistency");
          }
          for(unsigned int rr = 0;rr != rank;rr++) {
            cblas_daxpy(
              nb_row_dofs,w*(fun_val[lhs])[rr],&data.getN()(gg,0),1,&(Nf[lhs])[rr],rank
            );
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

  /** \brief Set operators
  */
  template<typename FUNEVAL>
  PetscErrorCode setOperatorsVolume(
    const std::string &field_name,
    Mat A,
    std::vector<Vec> &vec_F,
    FUNEVAL &function_evaluator
  ) {
    PetscFunctionBegin;
    //add operator to calculate F vector
    feVolume.getOpPtrVector().push_back(new OpApproxVolume<FUNEVAL>(field_name,A,vec_F,function_evaluator));
    //add operator to calculate A matrix
    // if(A) {
    //   feVolume.getOpPtrVector().push_back(new OpApproxVolume<FUNEVAL>(field_name,A,vec_F,function_evaluator));
    // }
    PetscFunctionReturn(0);
  }

  /** \brief Set operators
  */
  template<typename FUNEVAL>
  PetscErrorCode setOperatorsFace(
    const std::string &field_name,
    Mat A,
    std::vector<Vec> &vec_F,
    FUNEVAL &function_evaluator
  ) {
    PetscFunctionBegin;
    //add operator to calculate F vector
    feFace.getOpPtrVector().push_back(new OpApproxFace<FUNEVAL>(field_name,A,vec_F,function_evaluator));
    //add operator to calculate A matrix
    // if(A) {
    //   feFace.getOpPtrVector().push_back(new OpApproxFace<FUNEVAL>(field_name,A,vec_F,function_evaluator));
    // }
    PetscFunctionReturn(0);
  }

  /** \brief assemble matrix and vector
    */
  template<typename FUNEVAL>
  PetscErrorCode loopMatrixAndVectorVolume(
    const std::string &problem_name,
    const std::string &fe_name,
    const std::string &field_name,
    Mat A,std::vector<Vec> &vec_F,
    FUNEVAL &function_evaluator
  ) {
    PetscFunctionBegin;
    PetscErrorCode ierr;
    ierr = setOperatorsVolume(field_name,A,vec_F,function_evaluator); CHKERRQ(ierr);
    if(A) {
      ierr = MatZeroEntries(A); CHKERRQ(ierr);
    }
    //calculate and assemble
    ierr = mField.loop_finite_elements(problem_name,fe_name,feVolume);  CHKERRQ(ierr);
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

  template<typename FUNEVAL>
  DEPRECATED PetscErrorCode loopMatrixAndVector(
    const std::string &problem_name,
    const std::string &fe_name,
    const std::string &field_name,
    Mat A,std::vector<Vec> &vec_F,
    FUNEVAL &function_evaluator
  ) {
    PetscFunctionBegin;
    PetscErrorCode ierr;
    ierr = loopMatrixAndVectorVolume(
      problem_name,fe_name,field_name,A,vec_F,function_evaluator
    ); CHKERRQ(ierr);
    PetscFunctionReturn(0);
  }

};

#endif //__FILEDAPPROXIMATION_HPP
