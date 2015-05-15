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
    
    
    struct CommonData {
      int rank;
      ublas::matrix<double> D_mat;
    };
    CommonData commonData;
    
    struct CommonFunctions {
      PetscErrorCode shapeMat(int rank, unsigned int gg, DataForcesAndSurcesCore::EntData &col_data, ublas::matrix<FieldData> &N_mat) {
        PetscFunctionBegin;
        int shape_size=col_data.getN().size2();
        //      cout<<"shape_size = "<<shape_size<<endl;
        //      cout<<"rank = "<<rank<<endl;
        N_mat.resize(rank,shape_size*rank);    N_mat.clear();
        int gg1=0;       int kk=0;
        for(int ii=0; ii<shape_size; ii++){ //number of shape funcitons
          for(int jj=0; jj<rank; jj++){
            //              cout<<"ii "<<ii<<endl;
            N_mat(jj,kk)=col_data.getN()(gg,gg1); kk++;
          }
          gg1++;
        }
        PetscFunctionReturn(0);
      }
      
      
      PetscErrorCode HMat(ublas::vector<int> getNormals_at_GaussPt, int rank, ublas::matrix<FieldData> &N_mat, ublas::matrix<FieldData> &H_mat) {
        PetscFunctionBegin;
        //        cout<<"getNormals_at_GaussPt(gg) "<<getNormals_at_GaussPt(gg)<<endl;
        ublas::matrix<FieldData> H_mat_1Node;
        switch(rank) { //Mechanical problem
          case 3:
          {
            H_mat_1Node.resize(6,3);  H_mat_1Node.clear(); //for one node
            if(getNormals_at_GaussPt(0)>0){     //+X face of the RVE
              H_mat_1Node(0,0)=1.0;  H_mat_1Node(3,1)=1.0;  H_mat_1Node(4,2)=1.0;
            }
            if(getNormals_at_GaussPt(0)<0){    //-X face of the RVE
              H_mat_1Node(0,0)=-1.0;  H_mat_1Node(3,1)=-1.0;  H_mat_1Node(4,2)=-1.0;
            }
            if(getNormals_at_GaussPt(1)>0){     //+Y face of the RVE
              H_mat_1Node(1,1)=1.0;  H_mat_1Node(3,0)=1.0;  H_mat_1Node(5,2)=1.0;
            }
            if(getNormals_at_GaussPt(1)<0){    //-Y face of the RVE
              H_mat_1Node(1,1)=-1.0;  H_mat_1Node(3,0)=-1.0;  H_mat_1Node(5,2)=-1.0;
            }
            if(getNormals_at_GaussPt(2)>0){    //+Z face of the RVE
              H_mat_1Node(2,2)=1.0;  H_mat_1Node(4,0)=1.0;  H_mat_1Node(5,1)=1.0;
            }
            if(getNormals_at_GaussPt(2)<0){    //-Z face of the RVE
              H_mat_1Node(2,2)=-1.0;  H_mat_1Node(4,0)=-1.0;  H_mat_1Node(5,1)=-1.0;
            }
            
            int num_col=N_mat.size2(); H_mat.resize(6,num_col);
            int cc1=0;
            for(int bb = 0; bb<num_col/3; bb++) {  //blocks of 6x3
              for(int rr = 0; rr<6; rr++) {
                for(int cc = 0; cc<3; cc++) {
                  H_mat(rr,(cc+cc1))=H_mat_1Node(rr,cc);
                }
              }
              cc1+=3;
            }
            //            cout<<"H_mat "<<H_mat<<endl;
          }
            break;
          case 1: //Moisture transport or thermal problem
          {
            H_mat_1Node.resize(3,1);  H_mat_1Node.clear();
            if(getNormals_at_GaussPt(0)>0){     //+X face of the RVE
              H_mat_1Node(0,0)=1.0;
            }
            if(getNormals_at_GaussPt(0)<0){    //-X face of the RVE
              H_mat_1Node(0,0)=-1.0;
            }
            if(getNormals_at_GaussPt(1)>0){     //+Y face of the RVE
              H_mat_1Node(1,0)=1.0;
            }
            if(getNormals_at_GaussPt(1)<0){    //-Y face of the RVE
              H_mat_1Node(1,0)=-1.0;
            }
            if(getNormals_at_GaussPt(2)>0){    //+Z face of the RVE
              H_mat_1Node(2,0)=1.0;
            }
            if(getNormals_at_GaussPt(2)<0){    //-Z face of the RVE
              H_mat_1Node(2,0)=-1.0;
            }
            int num_col=N_mat.size2(); H_mat.resize(6,num_col);
            H_mat.resize(3,num_col);
            int cc1=0;
            for(int bb = 0; bb<num_col; bb++) {  //blocks of 3x1
              for(int rr = 0; rr<3; rr++) {
                for(int cc = 0; cc<1; cc++) {
                  H_mat(rr,(cc+cc1))=H_mat_1Node(rr,cc);
                }
              }
              cc1+=1;
            }
          }
            break;
          default:
            SETERRQ(PETSC_COMM_SELF,1,"not implemented");
        }
        PetscFunctionReturn(0);
      }
    };
    CommonFunctions common_functions;

    
    /// \biref operator to calculate the LHS for the RVE bounary conditions
    struct OpRVEBCsLhs:public FaceElementForcesAndSourcesCore::UserDataOperator  {
      
//      field_name, lagrang_field_name,
      CommonFunctions common_functions;
      RVEBC_Data &dAta;
      bool ho_geometry;
      Mat Aij;
      OpRVEBCsLhs(const string field_name, const string lagrang_field_name, Mat _Aij,
                  RVEBC_Data &data, CommonFunctions &_common_functions, bool _ho_geometry = false):
      FaceElementForcesAndSourcesCore::UserDataOperator(lagrang_field_name, field_name),Aij(_Aij),
      dAta(data),ho_geometry(_ho_geometry), common_functions(_common_functions){
        sYmm = false;  //This will make sure to loop over all intities (e.g. for order=2 it will make doWork to loop 16 time)
      }
      
      ublas::matrix<double> K,transK, NTN;
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
          
//          cout<<"OpRVEBCsLhs "<<endl;
          int nb_row = row_data.getIndices().size();
          int nb_col = col_data.getIndices().size();
          //          cout<<"nb_row "<< nb_row << endl;
          //          cout<<"nb_col "<< nb_col << endl;
          PetscErrorCode ierr;
          const FENumeredDofMoFEMEntity *dof_ptr;
          ierr = getMoFEMFEPtr()->get_col_dofs_by_petsc_gloabl_dof_idx(col_data.getIndices()[0],&dof_ptr); CHKERRQ(ierr);
          int rank = dof_ptr->get_max_rank();
          //          cout<<"rank "<< rank<< endl;
          
          K.resize(nb_row,nb_col);
          K.clear();
          //        cout<<"no of Gauss points per element = "<<row_data.getN().size1()<<endl;
          ublas::matrix<FieldData> N_mat, H_mat;
          for(unsigned int gg = 0;gg<col_data.getN().size1();gg++) {
            double area;
            if(ho_geometry) {
              area = norm_2(getNormals_at_GaussPt(gg))*0.5;
              //              cout<<"area with ho_geometry = "<<area<<endl;
            }   else {
              area = getArea();
              //              cout<<"area without ho_geometry = "<<area<<endl;
            }
            double val = getGaussPts()(2,gg)*area;
            
            ierr = common_functions.shapeMat(rank, gg, col_data, N_mat); CHKERRQ(ierr);
//            cout<<"N_mat "<<N_mat<<endl;
            
            ierr = common_functions.HMat(getNormals_at_GaussPt(gg), rank, N_mat, H_mat); CHKERRQ(ierr);
//            cout<<"H_mat "<<H_mat<<endl;
            if(gg==0){
              NTN=prod(trans(N_mat),N_mat);  //we don't need to define its size NTN
              K=val*prod(H_mat, NTN); //we don't need to defien size of K
            }else{
              NTN=prod(trans(N_mat),N_mat);
              K+=val*prod(H_mat, NTN);
            }
          }
//          cout<<"K = "<<K<<endl;
          // Matrix C
          int nb_rows=row_data.getIndices().size();
          int nb_cols=col_data.getIndices().size();
          ierr = MatSetValues(Aij,nb_rows,&row_data.getIndices()[0],nb_cols,&col_data.getIndices()[0],&K(0,0),ADD_VALUES); CHKERRQ(ierr);
          
          // Matrix CT
          transK.resize(N_mat.size2(), H_mat.size1());
          noalias(transK) = trans(K);
          ierr = MatSetValues(Aij,nb_cols,&col_data.getIndices()[0],nb_rows,&row_data.getIndices()[0],&transK(0,0),ADD_VALUES); CHKERRQ(ierr);
          
        } catch (const std::exception& ex) {
          ostringstream ss;
          ss << "throw in method: " << ex.what() << endl;
          SETERRQ(PETSC_COMM_SELF,1,ss.str().c_str());
        }
        
        PetscFunctionReturn(0);
      }
      
      
    };
    
    
    /// \biref operator to calculate the RHS of the constrain for the RVE boundary conditions
    struct OpRVEBCsRhs_Cal:public FaceElementForcesAndSourcesCore::UserDataOperator {
      
      RVEBC_Data &dAta;
      bool ho_geometry;
      CommonData &commonData;
      CommonFunctions common_functions;
      
      OpRVEBCsRhs_Cal(const string field_name, RVEBC_Data &data, CommonData &_commonData, CommonFunctions &_common_functions, bool _ho_geometry = false): FaceElementForcesAndSourcesCore::UserDataOperator(field_name),dAta(data), commonData(_commonData), common_functions(_common_functions),  ho_geometry(_ho_geometry){ }
      
      ublas::vector<FieldData> f;
      
      /** brief calculate Convection condition on the right hand side
       *  R=int_S N^T*alpha*N_f  dS **/
      PetscErrorCode doWork(int side,EntityType type,DataForcesAndSurcesCore::EntData &data) {
        PetscFunctionBegin;
        if(data.getIndices().size()==0) PetscFunctionReturn(0);
//        cout<<"data.getIndices().size() "<<data.getIndices().size()<<endl;
        PetscErrorCode ierr;
        const FENumeredDofMoFEMEntity *dof_ptr;
        ierr = getMoFEMFEPtr()->get_col_dofs_by_petsc_gloabl_dof_idx(data.getIndices()[0],&dof_ptr); CHKERRQ(ierr);
        int rank = dof_ptr->get_max_rank();
        commonData.rank=rank;
//        cout<<"commonData.rank "<<commonData.rank<<endl;

        ublas::matrix<FieldData> X_mat, N_mat, NTX, H_mat;
        X_mat.resize(rank,1.5*rank+1.5);    X_mat.clear();
        
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
              X_mat(0,0)=2.0*x;  X_mat(0,3)=y;  X_mat(0,4)=z;
              X_mat(1,1)=2.0*y;  X_mat(1,3)=x;  X_mat(1,5)=z;
              X_mat(2,2)=2.0*z;  X_mat(2,4)=x;  X_mat(2,5)=y;
              X_mat=0.5*X_mat;
              break;
            case 1:  //moisture transport or thermal problem
              X_mat(0,0)=x;  X_mat(0,1)=y;  X_mat(0,2)=z;
              break;
            default:
              SETERRQ(PETSC_COMM_SELF,1,"not implemented");
          }
          ierr = common_functions.shapeMat(rank, gg,  data, N_mat); CHKERRQ(ierr);
//          cout<<"N_mat "<<N_mat<<endl;
          ierr = common_functions.HMat(getNormals_at_GaussPt(gg), rank, N_mat, H_mat); CHKERRQ(ierr);
//          cout<<"H_mat "<<H_mat<<endl;

          if(gg==0){
            NTX=val*prod(trans(N_mat),X_mat);
            commonData.D_mat=prod(H_mat, NTX);
          }else{
            NTX=val*prod(trans(N_mat),X_mat);
            commonData.D_mat+=prod(H_mat, NTX);
          }
        }
//       cout<<"commonData.D_mat "<<commonData.D_mat<<endl;
       PetscFunctionReturn(0);
      }
    };

    
    
    
    
    /// \biref operator to calculate the RHS of the constrain for the RVE boundary conditions
    struct OpRVEBCsRhs_Assemble: public FaceElementForcesAndSourcesCore::UserDataOperator {
      
      RVEBC_Data &dAta;
      bool ho_geometry;
      CommonData &commonData;
      Vec F1, F2, F3, F4, F5, F6;
      
      OpRVEBCsRhs_Assemble(const string lagrang_field_name, Vec _F1, Vec _F2, Vec _F3, Vec _F4, Vec _F5, Vec _F6, RVEBC_Data &data, CommonData &_commonData, bool _ho_geometry = false): FaceElementForcesAndSourcesCore::UserDataOperator(lagrang_field_name),dAta(data), F1(_F1), F2(_F2), F3(_F3), F4(_F4), F5(_F5), F6(_F6), commonData(_commonData), ho_geometry(_ho_geometry){ }
      
      ublas::vector<FieldData> f;
      
      /** brief calculate Convection condition on the right hand side
       *  R=int_S N^T*alpha*N_f  dS **/
      PetscErrorCode doWork(int side,EntityType type,DataForcesAndSurcesCore::EntData &data) {
        PetscFunctionBegin;
//        cout<<"rank "<<commonData.rank<<endl;
//        cout<<"D_mat "<<commonData.D_mat<<endl;
        int rank =commonData.rank;
        ublas::matrix<FieldData> D_mat=commonData.D_mat;
//        ublas::vector<FieldData> applied_strain;
//        applied_strain.resize(1.5*rank+1.5);
//
////        ublas::vector<FieldData> applied_strain;
////        applied_strain.resize(1.5*rank+1.5);
//        applied_strain.clear();
//        applied_strain(0)=1.0;
//        //        cout<<"applied_strain = "<<applied_strain<<endl;
//        PetscErrorCode ierr;
//        f=prod(D_mat, applied_strain);
//        ierr = VecSetValues(F1,data.getIndices().size(),&data.getIndices()[0],&f[0],ADD_VALUES); CHKERRQ(ierr);
//        //      cout<<"Hello "<<endl;
//        //      std::string wait;
//        //      std::cin >> wait;
//        applied_strain.clear();
//        applied_strain(1)=1.0;
//        //      cout<<"applied_strain = "<<applied_strain<<endl;
//        f=prod(D_mat, applied_strain);
//        //      cout<<"f = "<<f<<endl;
//        ierr = VecSetValues(F2,data.getIndices().size(),&data.getIndices()[0],&f[0],ADD_VALUES); CHKERRQ(ierr);
//        
//        applied_strain.clear();
//        applied_strain(2)=1.0;
//        //      cout<<"applied_strain = "<<applied_strain<<endl;
//        f=prod(D_mat, applied_strain);
//        //      cout<<"f = "<<f<<endl;
//        ierr = VecSetValues(F3,data.getIndices().size(),&data.getIndices()[0],&f[0],ADD_VALUES); CHKERRQ(ierr);
//        
//        if(rank==3){ //This poriton will execute for mechanical problem (rank=3) only
//          applied_strain.clear();
//          applied_strain(3)=1.0;
//          //      cout<<"applied_strain = "<<applied_strain<<endl;
//          f=prod(D_mat, applied_strain);
//          //      cout<<"f = "<<f<<endl;
//          ierr = VecSetValues(F4,data.getIndices().size(),&data.getIndices()[0],&f[0],ADD_VALUES); CHKERRQ(ierr);
//          
//          applied_strain.clear();
//          applied_strain(4)=1.0;
//          //      cout<<"applied_strain = "<<applied_strain<<endl;
//          f=prod(D_mat, applied_strain);
//          //      cout<<"f = "<<f<<endl;
//          ierr = VecSetValues(F5,data.getIndices().size(),&data.getIndices()[0],&f[0],ADD_VALUES); CHKERRQ(ierr);
//          
//          applied_strain.clear();
//          applied_strain(5)=1.0;
//          //      cout<<"applied_strain = "<<applied_strain<<endl;
//          f=prod(D_mat, applied_strain);
//          //      cout<<"f = "<<f<<endl;
//          ierr = VecSetValues(F6,data.getIndices().size(),&data.getIndices()[0],&f[0],ADD_VALUES); CHKERRQ(ierr);
//        }
        PetscFunctionReturn(0);
      }
    };

    
    
     PetscErrorCode setRVEBCsOperators(string field_name,string lagrang_field_name,Mat _Aij, Vec _F1, Vec _F2, Vec _F3, Vec _F4, Vec _F5, Vec _F6, const string mesh_nodals_positions) {
      PetscFunctionBegin;
      
      bool ho_geometry = false;
      if(mField.check_field(mesh_nodals_positions)) {
        ho_geometry = true;
      }
//      cout<<"Hi 1 from setRVEBCsOperators "<<endl;
      map<int,RVEBC_Data>::iterator sit = setOfRVEBC.begin();
      for(;sit!=setOfRVEBC.end();sit++) {
//        cout<<"Hi from setRVEBCsOperators "<<endl;
        //LHS
        feRVEBCLhs.getRowColOpPtrVector().push_back(new OpRVEBCsLhs(field_name,lagrang_field_name, _Aij, sit->second, common_functions, ho_geometry));
        
        //RHS
        
//        feRVEBCRhs.getOpPtrVector().push_back(new OpRVEBCsRhs_Cal(field_name, sit->second, commonData, common_functions, ho_geometry,ForcesAndSurcesCore::UserDataOperator::OPROW));
//        fe1.getOpPtrVector().push_back(new MyOp1(field_name_col,ForcesAndSurcesCore::UserDataOperator::OPCOL));
        
        
        
        //Caclculte D_mat
        feRVEBCRhs.getColOpPtrVector().push_back(new OpRVEBCsRhs_Cal(field_name, sit->second, commonData, common_functions, ho_geometry));
        //Caclculte f and assemplbe
//        feRVEBCRhs.getRowOpPtrVector().push_back(new OpRVEBCsRhs_Assemble(lagrang_field_name, _F1, _F2, _F3, _F4, _F5, _F6, sit->second, commonData, ho_geometry));
      }
      PetscFunctionReturn(0);
    }
    
    
    
  };
  
}

#endif
