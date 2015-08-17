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

#ifndef __BCS_RVELAGRANGE_DISP_HPP
#define __BCS_RVELAGRANGE_DISP_HPP

namespace MoFEM {
  
  struct BCs_RVELagrange_Disp {
    
    FieldInterface &mField;
    
    struct MyTriFE: public FaceElementForcesAndSourcesCore {
      MyTriFE(FieldInterface &_mField): FaceElementForcesAndSourcesCore(_mField) {}
      int getRule(int order) { return order; };
    };

    MyTriFE feRVEBCLhs; //To calculate the Lhs or RVE BCs
    MyTriFE feRVEBCRhs; //To calculate the Rhs or RVE BCs
    
    MyTriFE& getLoopFeRVEBCLhs() { return feRVEBCLhs; }
    MyTriFE& getLoopFeRVEBCRhs() { return feRVEBCRhs; }


    BCs_RVELagrange_Disp(FieldInterface &m_field):
    feRVEBCRhs(m_field),feRVEBCLhs(m_field),
    mField(m_field) {}

    struct RVEBC_Data {
      Range tRis; // All boundary surfaces
    };
    map<int,RVEBC_Data> setOfRVEBC; ///< maps side set id with appropriate FluxData
    
    
    PetscErrorCode addLagrangiangElement(const string element_name, const string field_name, const string lagrang_field_name, const string mesh_nodals_positions) {
      PetscFunctionBegin;
      
      PetscErrorCode ierr;
      ErrorCode rval;
      
      ierr = mField.add_finite_element(element_name,MF_ZERO); CHKERRQ(ierr);
      //============================================================================================================
      //C row as Lagrange_mul_disp and col as DISPLACEMENT
      ierr = mField.modify_finite_element_add_field_row(element_name,lagrang_field_name); CHKERRQ(ierr);
      ierr = mField.modify_finite_element_add_field_col(element_name,field_name); CHKERRQ(ierr);
      //CT col as Lagrange_mul_disp and row as DISPLACEMENT
      ierr = mField.modify_finite_element_add_field_col(element_name,lagrang_field_name); CHKERRQ(ierr);
      ierr = mField.modify_finite_element_add_field_row(element_name,field_name); CHKERRQ(ierr);
      //data
      ierr = mField.modify_finite_element_add_field_data(element_name,lagrang_field_name); CHKERRQ(ierr);
      ierr = mField.modify_finite_element_add_field_data(element_name,field_name); CHKERRQ(ierr);
      //============================================================================================================

      if(mField.check_field(mesh_nodals_positions)) { //for high order geometry
        ierr = mField.modify_finite_element_add_field_data(element_name,mesh_nodals_positions); CHKERRQ(ierr);
      }
      //this is alternative method for setting boundary conditions, to bypass bu in cubit file reader.
      //not elegant, but good enough
      for(_IT_CUBITMESHSETS_BY_SET_TYPE_FOR_LOOP_(mField,SIDESET,it)) {
        if(it->get_name().compare(0,12,"AllBoundSurf") == 0) {
//          cout<<"Hi from addLagrangiangElement "<<endl;
          rval = mField.get_moab().get_entities_by_type(it->meshset,MBTRI,setOfRVEBC[it->get_msId()].tRis,true); CHKERR_PETSC(rval);
          ierr = mField.add_ents_to_finite_element_by_TRIs(setOfRVEBC[it->get_msId()].tRis,element_name); CHKERRQ(ierr);
          
        }
      }
       PetscFunctionReturn(0);
    }

    
    
    
    /// \biref operator to calculate the LHS for the RVE bounary conditions
    struct OpRVEBCsLhs:public FaceElementForcesAndSourcesCore::UserDataOperator {
      
      RVEBC_Data &dAta;
      bool ho_geometry;
      Mat Aij;
      OpRVEBCsLhs(const string field_name, const string lagrang_field_name, Mat _Aij,
                      RVEBC_Data &data,bool _ho_geometry = false):
      FaceElementForcesAndSourcesCore::UserDataOperator(lagrang_field_name, field_name, UserDataOperator::OPROWCOL),Aij(_Aij),
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
//          cout<<"OpRVEBCsLhs "<<endl;
          if(row_data.getIndices().size()==0) PetscFunctionReturn(0);
          if(col_data.getIndices().size()==0) PetscFunctionReturn(0);
          
          int nb_row = row_data.getIndices().size();
          int nb_col = col_data.getIndices().size();
//          cout<<"nb_row "<< nb_row << endl;
//          cout<<"nb_col "<< nb_col << endl;
          PetscErrorCode ierr;
          const FENumeredDofMoFEMEntity *dof_ptr;
          ierr = getMoFEMFEPtr()->get_row_dofs_by_petsc_gloabl_dof_idx(row_data.getIndices()[0],&dof_ptr); CHKERRQ(ierr);
          int rank = dof_ptr->get_max_rank();
//          cout<<"rank "<< rank<< endl;

          K.resize(nb_row/rank,nb_col/rank);
          K.clear();
//          cout<<"no of Gauss points per element = "<<row_data.getN().size1()<<endl;

          for(unsigned int gg = 0;gg<row_data.getN().size1();gg++) {
            double area;
            if(ho_geometry) {
              area = norm_2(getNormals_at_GaussPt(gg))*0.5;
//              cout<<"area with ho_geometry = "<<area<<endl;
            }   else {
              area = getArea();
//              cout<<"area without ho_geometry = "<<area<<endl;
            }
            double val = getGaussPts()(2,gg)*area;
            noalias(K) += val*outer_prod(row_data.getN(gg,nb_row/rank),col_data.getN(gg,nb_col/rank) );
          }
//          cout<<"K "<< K<< endl;
//          cout<<"&row_data.getIndices()[0] "<< endl;
//          for (int ii=0; ii<nb_row; ii++) cout<<row_data.getIndices()[ii]<<"  ";
//          cout<< endl;
//          cout<<"&col_data.getIndices()[0] "<< endl;
//          for (int ii=0; ii<nb_col; ii++) cout<<col_data.getIndices()[ii]<<"  ";
//          cout<< endl<<endl;
          ublas::vector<DofIdx> row_indices,col_indices;
          row_indices.resize(nb_row/rank);
          col_indices.resize(nb_col/rank);

          for(int rr = 0;rr < rank; rr++) {
            unsigned int nb_rows;
            unsigned int nb_cols;
            int *rows;
            int *cols;
            if(rank > 1) {
              ublas::noalias(row_indices) = ublas::vector_slice<ublas::vector<DofIdx> >
              (row_data.getIndices(), ublas::slice(rr, rank, row_data.getIndices().size()/rank));
//              cout<<"&row_indices "<< row_indices << endl;

              ublas::noalias(col_indices) = ublas::vector_slice<ublas::vector<DofIdx> >
              (col_data.getIndices(), ublas::slice(rr, rank, col_data.getIndices().size()/rank));
//              cout<<"&col_indices "<< col_indices << endl;

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

            // Matrix C
            ierr = MatSetValues(Aij,nb_rows,rows,nb_cols,cols,&K(0,0),ADD_VALUES); CHKERRQ(ierr);
            
            // Matrix CT
            transK.resize(nb_col/rank,nb_row/rank);
            noalias(transK) = trans(K);
            ierr = MatSetValues(Aij,nb_cols,cols,nb_rows,rows,&transK(0,0),ADD_VALUES); CHKERRQ(ierr);
          }
          
        } catch (const std::exception& ex) {
          ostringstream ss;
          ss << "throw in method: " << ex.what() << endl;
          SETERRQ(PETSC_COMM_SELF,1,ss.str().c_str());
        }
        
        PetscFunctionReturn(0);
      }
    };

    
    
  
    
    /// \biref operator to calculate the RHS of the constrain for the RVE boundary conditions
    struct OpRVEBCsRhs:public FaceElementForcesAndSourcesCore::UserDataOperator {
      
      RVEBC_Data &dAta;
      bool ho_geometry;
      Vec F1, F2, F3, F4, F5, F6;
      
      OpRVEBCsRhs(const string field_name, const string lagrang_field_name,  Vec _F1, Vec _F2, Vec _F3, Vec _F4, Vec _F5, Vec _F6,
                  RVEBC_Data &data, bool _ho_geometry = false):
      FaceElementForcesAndSourcesCore::UserDataOperator(lagrang_field_name, UserDataOperator::OPROW),F1(_F1), F2(_F2), F3(_F3), F4(_F4), F5(_F5), F6(_F6),
      dAta(data),ho_geometry(_ho_geometry){}
      
      ublas::vector<FieldData> f;
      
      /** brief calculate Convection condition on the right hand side
       *  R=int_S N^T*alpha*N_f  dS **/
      PetscErrorCode doWork(int side,EntityType type,DataForcesAndSurcesCore::EntData &data) {
        PetscFunctionBegin;
        if(data.getIndices().size()==0) PetscFunctionReturn(0);
        if(dAta.tRis.find(getMoFEMFEPtr()->get_ent())==dAta.tRis.end()) PetscFunctionReturn(0);
//        cout<<"OpRVEBCsRhs "<<endl;

        PetscErrorCode ierr;
        const FENumeredDofMoFEMEntity *dof_ptr;
        ierr = getMoFEMFEPtr()->get_row_dofs_by_petsc_gloabl_dof_idx(data.getIndices()[0],&dof_ptr); CHKERRQ(ierr);
        int rank = dof_ptr->get_max_rank();
        int nb_row_dofs = data.getIndices().size()/rank;
        
        ublas::vector<FieldData> f;

        ublas::matrix<FieldData> X_mat, N_mat;
        X_mat.resize(rank,1.5*rank+1.5);    X_mat.clear();
        ublas::matrix<FieldData>  D_mat;
        
        ublas::vector<FieldData> applied_strain;
        applied_strain.resize(1.5*rank+1.5);


        double x,y,z;
        for(unsigned int gg = 0;gg<data.getN().size1();gg++) {
          double area;
          if(ho_geometry) {
            area = norm_2(getNormals_at_GaussPt(gg))*0.5;
          } else {
            area = getArea();
          }
          double val = getGaussPts()(2,gg)*area;

          x = getCoordsAtGaussPts()(gg,0);
          y = getCoordsAtGaussPts()(gg,1);
          z = getCoordsAtGaussPts()(gg,2);
          //          cout<<"Gauss point coordinates "<< "  x=    "<< x << "   y=   "<< y << "  z=    "<<z<<endl;

          switch(rank) {
            case 3: //mech problem
              X_mat(0,0)=2.0*x;  X_mat(0,3)=y;  X_mat(0,5)=z;
              X_mat(1,1)=2.0*y;  X_mat(1,3)=x;  X_mat(1,4)=z;
              X_mat(2,2)=2.0*z;  X_mat(2,4)=y;  X_mat(2,5)=x;
              X_mat=0.5*X_mat;
              break;
            case 1:  //moisture transport or thermal problem
              X_mat(0,0)=x;  X_mat(0,1)=y;  X_mat(0,2)=z;
              break;
            default:
              SETERRQ(PETSC_COMM_SELF,1,"not implemented");
          }
          
          int shape_size=data.getN().size2();
          N_mat.resize(rank,shape_size*rank);   N_mat.clear();
          int gg1=0;       int kk=0;
          for(int ii=0; ii<shape_size; ii++){ //number of shape funcitons
            for(int jj=0; jj<rank; jj++){
//              cout<<"ii "<<ii<<endl;
              N_mat(jj,kk)=data.getN()(gg,gg1); kk++;
            }
            gg1++;
          }
//          cout<<"N_mat "<<N_mat<<endl;
          if(gg==0){
            D_mat=val*prod(trans(N_mat),X_mat);
          }else{
            D_mat+=val*prod(trans(N_mat),X_mat);
          }
        }
        
        
        applied_strain.clear();
        applied_strain(0)=1.0;
//        cout<<"applied_strain = "<<applied_strain<<endl;
        f=prod(D_mat, applied_strain);
        ierr = VecSetValues(F1,data.getIndices().size(),&data.getIndices()[0],&f[0],ADD_VALUES); CHKERRQ(ierr);
        //      cout<<"Hello "<<endl;
        //      std::string wait;
        //      std::cin >> wait;
        applied_strain.clear();
        applied_strain(1)=1.0;
        //      cout<<"applied_strain = "<<applied_strain<<endl;
        f=prod(D_mat, applied_strain);
        //      cout<<"f = "<<f<<endl;
        ierr = VecSetValues(F2,data.getIndices().size(),&data.getIndices()[0],&f[0],ADD_VALUES); CHKERRQ(ierr);
        
        applied_strain.clear();
        applied_strain(2)=1.0;
        //      cout<<"applied_strain = "<<applied_strain<<endl;
        f=prod(D_mat, applied_strain);
        //      cout<<"f = "<<f<<endl;
        ierr = VecSetValues(F3,data.getIndices().size(),&data.getIndices()[0],&f[0],ADD_VALUES); CHKERRQ(ierr);
        
        if(rank==3){ //This poriton will execute for mechanical problem (rank=3) only
          applied_strain.clear();
          applied_strain(3)=1.0;
          //      cout<<"applied_strain = "<<applied_strain<<endl;
          f=prod(D_mat, applied_strain);
          //      cout<<"f = "<<f<<endl;
          ierr = VecSetValues(F4,data.getIndices().size(),&data.getIndices()[0],&f[0],ADD_VALUES); CHKERRQ(ierr);
          
          applied_strain.clear();
          applied_strain(4)=1.0;
          //      cout<<"applied_strain = "<<applied_strain<<endl;
          f=prod(D_mat, applied_strain);
          //      cout<<"f = "<<f<<endl;
          ierr = VecSetValues(F5,data.getIndices().size(),&data.getIndices()[0],&f[0],ADD_VALUES); CHKERRQ(ierr);
          
          applied_strain.clear();
          applied_strain(5)=1.0;
          //      cout<<"applied_strain = "<<applied_strain<<endl;
          f=prod(D_mat, applied_strain);
          //      cout<<"f = "<<f<<endl;
          ierr = VecSetValues(F6,data.getIndices().size(),&data.getIndices()[0],&f[0],ADD_VALUES); CHKERRQ(ierr);
        }
//        f=prod(D_mat, applied_strain);
////        cout<<"f "<<f<<endl;
//        ierr = VecSetValues(F,data.getIndices().size(), &data.getIndices()[0],&f[0],ADD_VALUES); CHKERRQ(ierr);
        PetscFunctionReturn(0);
      }
     };

    
     
    PetscErrorCode setRVEBCsOperators(string field_name,string lagrang_field_name,Mat _Aij, Vec _F1, Vec _F2, Vec _F3, Vec _F4, Vec _F5, Vec _F6, const string mesh_nodals_positions) {
      PetscFunctionBegin;
      
      bool ho_geometry = false;
      if(mField.check_field(mesh_nodals_positions)) {
//        cout<<"ho_geometry checking "<<endl;
        ho_geometry = true;
      }
      //      cout<<"Hi 1 from setRVEBCsOperators "<<endl;
      map<int,RVEBC_Data>::iterator sit = setOfRVEBC.begin();
      for(;sit!=setOfRVEBC.end();sit++) {
        //        cout<<"Hi from setRVEBCsOperators "<<endl;
        feRVEBCLhs.getOpPtrVector().push_back(new OpRVEBCsLhs(field_name,lagrang_field_name, _Aij, sit->second,ho_geometry));
        feRVEBCRhs.getOpPtrVector().push_back(new OpRVEBCsRhs(field_name,lagrang_field_name, _F1, _F2, _F3, _F4, _F5, _F6, sit->second,ho_geometry));
        
      }
      PetscFunctionReturn(0);
    }
    

    
  };

}

#endif
