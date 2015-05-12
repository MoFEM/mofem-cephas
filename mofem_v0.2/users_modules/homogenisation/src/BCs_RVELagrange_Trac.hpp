/* Copyright (C) 2014, Zahur Ullah (Zahur.Ullah AT glasgow.ac.uk)
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

#ifndef __BCS_RVELAGRANGE_TRAC_HPP
#define __BCS_RVELAGRANGE_TRAC_HPP

namespace MoFEM {
  
  struct BCs_RVELagrange_Trac: public BCs_RVELagrange_Disp {
    
    BCs_RVELagrange_Trac(FieldInterface &m_field): BCs_RVELagrange_Disp(m_field){}

    
    
    /// \biref operator to calculate the LHS for the RVE bounary conditions
    struct OpRVEBCsLhs:public FaceElementForcesAndSourcesCore::UserDataOperator {
      
      RVEBC_Data &dAta;
      bool ho_geometry;
      Mat Aij;
      OpRVEBCsLhs(const string field_name, const string lagrang_field_name, Mat _Aij,
                  RVEBC_Data &data,bool _ho_geometry = false):
      FaceElementForcesAndSourcesCore::UserDataOperator(lagrang_field_name, field_name),Aij(_Aij),
      dAta(data),ho_geometry(_ho_geometry){
        sYmm = false;  //This will make sure to loop over all intities (e.g. for order=2 it will make doWork to loop 16 time)
      }
      
      ublas::matrix<double> K,transK;
      /** \brief calculate thermal convection term in the lhs of equations
       *
       * C = intS N^T  N dS
       */
      PetscErrorCode doWork(
                            int row_side,int col_side,
                            EntityType row_type,EntityType col_type,
                            DataForcesAndSurcesCore::EntData &row_data,
                            DataForcesAndSurcesCore::EntData &col_data) {
        PetscFunctionBegin;
        
        try {
          cout<<"OpRVEBCsLhs "<<endl;
//          if(row_data.getIndices().size()==0) PetscFunctionReturn(0);
//          if(col_data.getIndices().size()==0) PetscFunctionReturn(0);
//          
//          int nb_row = row_data.getIndices().size();
//          int nb_col = col_data.getIndices().size();
//          //          cout<<"nb_row "<< nb_row << endl;
//          //          cout<<"nb_col "<< nb_col << endl;
//          PetscErrorCode ierr;
//          const FENumeredDofMoFEMEntity *dof_ptr;
//          ierr = getMoFEMFEPtr()->get_row_dofs_by_petsc_gloabl_dof_idx(row_data.getIndices()[0],&dof_ptr); CHKERRQ(ierr);
//          int rank = dof_ptr->get_max_rank();
//          //          cout<<"rank "<< rank<< endl;
//          
//          K.resize(nb_row/rank,nb_col/rank);
//          K.clear();
//          //          cout<<"no of Gauss points per element = "<<row_data.getN().size1()<<endl;
//          
//          for(unsigned int gg = 0;gg<row_data.getN().size1();gg++) {
//            double area;
//            if(ho_geometry) {
//              area = norm_2(getNormals_at_GaussPt(gg))*0.5;
//              //              cout<<"area with ho_geometry = "<<area<<endl;
//            }   else {
//              area = getArea();
//              //              cout<<"area without ho_geometry = "<<area<<endl;
//            }
//            
//            
//            double val = getGaussPts()(2,gg)*area;
//            noalias(K) += val*outer_prod(row_data.getN(gg,nb_row/rank),col_data.getN(gg,nb_col/rank) );
//          }
//          //          cout<<"K "<< K<< endl;
//          //          cout<<"&row_data.getIndices()[0] "<< endl;
//          //          for (int ii=0; ii<nb_row; ii++) cout<<row_data.getIndices()[ii]<<"  ";
//          //          cout<< endl;
//          //          cout<<"&col_data.getIndices()[0] "<< endl;
//          //          for (int ii=0; ii<nb_col; ii++) cout<<col_data.getIndices()[ii]<<"  ";
//          //          cout<< endl<<endl;
//          ublas::vector<DofIdx> row_indices,col_indices;
//          row_indices.resize(nb_row/rank);
//          col_indices.resize(nb_col/rank);
//          
//          for(int rr = 0;rr < rank; rr++) {
//            unsigned int nb_rows;
//            unsigned int nb_cols;
//            int *rows;
//            int *cols;
//            if(rank > 1) {
//              ublas::noalias(row_indices) = ublas::vector_slice<ublas::vector<DofIdx> >
//              (row_data.getIndices(), ublas::slice(rr, rank, row_data.getIndices().size()/rank));
//              //              cout<<"&row_indices "<< row_indices << endl;
//              
//              ublas::noalias(col_indices) = ublas::vector_slice<ublas::vector<DofIdx> >
//              (col_data.getIndices(), ublas::slice(rr, rank, col_data.getIndices().size()/rank));
//              //              cout<<"&col_indices "<< col_indices << endl;
//              
//              nb_rows = row_indices.size();
//              nb_cols = col_indices.size();
//              rows = &*row_indices.data().begin();
//              cols = &*col_indices.data().begin();
//            } else {
//              nb_rows = row_data.getIndices().size();
//              nb_cols = col_data.getIndices().size();
//              rows = &*row_data.getIndices().data().begin();
//              cols = &*col_data.getIndices().data().begin();
//            }
//            
//            // Matrix C
//            ierr = MatSetValues(Aij,nb_rows,rows,nb_cols,cols,&K(0,0),ADD_VALUES); CHKERRQ(ierr);
//            
//            // Matrix CT
//            transK.resize(nb_col/rank,nb_row/rank);
//            noalias(transK) = trans(K);
//            ierr = MatSetValues(Aij,nb_cols,cols,nb_rows,rows,&transK(0,0),ADD_VALUES); CHKERRQ(ierr);
//          }
//          
        } catch (const std::exception& ex) {
          ostringstream ss;
          ss << "throw in method: " << ex.what() << endl;
          SETERRQ(PETSC_COMM_SELF,1,ss.str().c_str());
        }
        
        PetscFunctionReturn(0);
      }
    };

    
    
    
    
    
    
    PetscErrorCode setRVEBCsOperators(string field_name,string lagrang_field_name,Mat _Aij, Vec _F1, Vec _F2, Vec _F3, Vec _F4, Vec _F5, Vec _F6, const string mesh_nodals_positions) {
      PetscFunctionBegin;
      
      bool ho_geometry = false;
      if(mField.check_field(mesh_nodals_positions)) {
        ho_geometry = true;
      }
      cout<<"Hi 1 from setRVEBCsOperators "<<endl;
      map<int,RVEBC_Data>::iterator sit = setOfRVEBC.begin();
      for(;sit!=setOfRVEBC.end();sit++) {
        cout<<"Hi from setRVEBCsOperators "<<endl;
        feRVEBCLhs.getRowColOpPtrVector().push_back(new OpRVEBCsLhs(field_name,lagrang_field_name, _Aij, sit->second,ho_geometry));
//        feRVEBCRhs.getRowOpPtrVector().push_back(new OpRVEBCsRhs(field_name,lagrang_field_name, _F1, _F2, _F3, _F4, _F5, _F6, sit->second,ho_geometry));
        
      }
      PetscFunctionReturn(0);
    }

    
    
  };

}

#endif
