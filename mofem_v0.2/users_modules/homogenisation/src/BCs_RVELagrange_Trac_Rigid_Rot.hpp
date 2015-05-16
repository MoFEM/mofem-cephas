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

#ifndef __BCS_RVELAGRANGE_TRAC_RIGID_ROT_HPP
#define __BCS_RVELAGRANGE_TRAC_RIGID_ROT_HPP

namespace MoFEM {
  
  struct BCs_RVELagrange_Trac_Rigid_Rot: public BCs_RVELagrange_Disp {
    
    BCs_RVELagrange_Trac_Rigid_Rot(FieldInterface &m_field): BCs_RVELagrange_Disp(m_field){}
    
    
    
    
    /// \biref operator to calculate the LHS for the RVE bounary conditions
    struct OpRVEBCsLhs:public FaceElementForcesAndSourcesCore::UserDataOperator  {
      RVEBC_Data &dAta;
      bool ho_geometry;
      Mat Aij;
      FieldInterface &mField;

      OpRVEBCsLhs(FieldInterface &_mField, const string field_name, const string lagrang_field_name, Mat _Aij, RVEBC_Data &data):
      FaceElementForcesAndSourcesCore::UserDataOperator(lagrang_field_name, field_name, UserDataOperator::OPROWCOL),Aij(_Aij), mField(_mField), dAta(data){
        sYmm = false;  //This will make sure to loop over all intities (e.g. for order=2 it will make doWork to loop 16 time)
      }
      
      //      ublas::matrix<double> K,transK, NTN;
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
          if(row_data.getIndices().size()==0) PetscFunctionReturn(0);
          if(col_data.getIndices().size()==0) PetscFunctionReturn(0);
//          cout<<"row_data.getIndices().size() "<<row_data.getIndices().size()<<endl;
//          cout<<"col_data.getIndices().size() "<<col_data.getIndices().size()<<endl;
          PetscErrorCode ierr;
          ublas::matrix<FieldData> Mat_face;          Mat_face.resize(3,9);           Mat_face.clear();
          ublas::matrix<FieldData> Mat_face_Tran;     Mat_face_Tran.resize(9,3);      Mat_face_Tran.clear();

          
          int num_nodes;   const EntityHandle* conn_face;  double coords_face[9];
          EntityHandle face_tri = getMoFEMFEPtr()->get_ent(); //handle of finite element
//          EntityHandle face_tri;  face_tri=fePtr->get_ent();
          ErrorCode rval;
          rval = mField.get_moab().get_connectivity(face_tri,conn_face,num_nodes,true); CHKERR_PETSC(rval);
          rval = mField.get_moab().get_coords(conn_face,num_nodes,coords_face); CHKERR_PETSC(rval);
          
//          for(int ii=0; ii<9; ii++) cout<<"coord "<<coords_face[ii]<<endl;
//          cout<< "num_nodes ="<<num_nodes << endl;
//          cout<< "conn_face ="<<conn_face << endl;
          for(int nn=0; nn<3; nn++){
            Mat_face(0,3*nn+1)=-coords_face[3*nn+2];    Mat_face(0,3*nn+2)= coords_face[3*nn+1];
            Mat_face(1,3*nn+0)= coords_face[3*nn+2];    Mat_face(1,3*nn+2)=-coords_face[3*nn+0];
            Mat_face(2,3*nn+0)=-coords_face[3*nn+1];    Mat_face(2,3*nn+1)= coords_face[3*nn+0];
          }

          // Matrix C1
          int nb_rows=row_data.getIndices().size();
          int nb_cols=col_data.getIndices().size();
          ierr = MatSetValues(Aij,nb_rows,&row_data.getIndices()[0],nb_cols,&col_data.getIndices()[0],&Mat_face(0,0),ADD_VALUES); CHKERRQ(ierr);

          // Matrix C1T
          noalias(Mat_face_Tran) = trans(Mat_face);
          ierr = MatSetValues(Aij,nb_cols,&col_data.getIndices()[0],nb_rows,&row_data.getIndices()[0],&Mat_face_Tran(0,0),ADD_VALUES); CHKERRQ(ierr);
          
        } catch (const std::exception& ex) {
          ostringstream ss;
          ss << "throw in method: " << ex.what() << endl;
          SETERRQ(PETSC_COMM_SELF,1,ss.str().c_str());
        }
        PetscFunctionReturn(0);
      }
    };
    
    
    
    PetscErrorCode setRVEBCsRigidBodyRotOperators(string field_name,string lagrang_field_name,Mat _Aij, map<int,RVEBC_Data> setOfRVEBC) {
      PetscFunctionBegin;
//      cout<<"Hi 1 from setRVEBCsRigidBodyRotOperators "<<endl;
      map<int,RVEBC_Data>::iterator sit = setOfRVEBC.begin();
      for(;sit!=setOfRVEBC.end();sit++) {
//        cout<<"Hi from setRVEBCsOperators "<<endl;
        //LHS
        feRVEBCLhs.getOpPtrVector().push_back(new OpRVEBCsLhs(mField, field_name,lagrang_field_name, _Aij, sit->second));
      }
      PetscFunctionReturn(0);
    }
    
    
    
  };
  
}

#endif
