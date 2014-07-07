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

#ifndef __MOISTUREFEMETHOD_HPP__
#define __MOISTUREFEMETHOD_HPP__

#include <boost/numeric/ublas/symmetric.hpp>

namespace MoFEM {
  
  struct MoistureFEMethod: public FEMethod_UpLevelStudent {
    
    FieldInterface& mField;
    Mat Aij;
    Vec Data,F;
    
    MoistureFEMethod(
                     FieldInterface& _mField): FEMethod_UpLevelStudent(_mField.get_moab(),1), mField(_mField),
    Aij(PETSC_NULL),Data(PETSC_NULL),F(PETSC_NULL) {};
    
    MoistureFEMethod(
                     FieldInterface& _mField,BaseDirihletBC *_dirihlet_ptr,Mat &_Aij,Vec &_D,Vec& _F):
    FEMethod_UpLevelStudent(_mField.get_moab(),_dirihlet_ptr,1), mField(_mField),
    Aij(_Aij),Data(_D),F(_F){
      
      pcomm = ParallelComm::get_pcomm(&moab,MYPCOMM_INDEX);
      
      RowGlob.resize(1+6+4+1);
      RowLocal.resize(1+6+4+1);
      rowNMatrices.resize(1+6+4+1);
      rowDiffNMatrices.resize(1+6+4+1);
      
      ColGlob.resize(1+6+4+1);
      ColLocal.resize(1+6+4+1);
      colNMatrices.resize(1+6+4+1);
      colDiffNMatrices.resize(1+6+4+1);
      
      g_NTET.resize(4*45);
      ShapeMBTET(&g_NTET[0],G_TET_X45,G_TET_Y45,G_TET_Z45,45);
      G_W_TET = G_TET_W45;
    };
    
    ErrorCode rval;
    ParallelComm* pcomm;
    PetscLogDouble t1,t2;
    PetscLogDouble v1,v2;
    
    
    
    int row_mat,col_mat;
    vector<vector<DofIdx> > RowGlob;
    vector<vector<DofIdx> > RowLocal;
    vector<vector<ublas::matrix<FieldData> > > rowNMatrices;
    vector<vector<ublas::matrix<FieldData> > > rowDiffNMatrices;
    
    vector<vector<DofIdx> > ColGlob;
    vector<vector<DofIdx> > ColLocal;
    vector<vector<ublas::matrix<FieldData> > > colNMatrices;
    vector<vector<ublas::matrix<FieldData> > > colDiffNMatrices;
    
    vector<double> g_NTET,g_NTRI;
    const double* G_W_TET;
    const double* G_W_TRI;
    
    
    ublas::matrix<FieldData> D;
    virtual PetscErrorCode calculateD(double _Moist_diffusivity) {
      PetscFunctionBegin;
      
      ublas::matrix<FieldData> D_lambda;
      D_lambda.resize(3,3);
      D_lambda.clear();
      for(int rr = 0;rr<3;rr++) {
        for(int cc = 0;cc<3;cc++) {
          if(rr==cc) D_lambda(rr,cc) = 1;
        }
      }
      
      //      cout<<"D_lambda "<<D_lambda<<endl;
      D = _Moist_diffusivity*D_lambda;
      //      cout<<"D "<<D<<endl;
      
      PetscFunctionReturn(0);
    }
    
    PetscErrorCode GetMatParameters(double *_Moist_diffusivity) {
      PetscFunctionBegin;
      
      EntityHandle ent = fe_ptr->get_ent();
      for(_IT_CUBITMESHSETS_BY_BCDATA_TYPE_FOR_LOOP_(mField,BlockSet|Mat_MoistureSet,it)) {
        
        Mat_Moisture mydata;
        ierr = it->get_attribute_data_structure(mydata); CHKERRQ(ierr);
        
        Range meshsets;
        rval = moab.get_entities_by_type(it->meshset,MBENTITYSET,meshsets,true); CHKERR_PETSC(rval);
        for(Range::iterator mit = meshsets.begin();mit != meshsets.end(); mit++) {
          if( moab.contains_entities(*mit,&ent,1) ) {
            *_Moist_diffusivity = mydata.data.Diffusivity;
            //                    cout<< " mydata.data.Diffusivity "<<mydata.data.Diffusivity<<endl;
            break;
          }
        }
      }
      PetscFunctionReturn(0);
    }
    
    
    
    PetscErrorCode preProcess() {
      PetscFunctionBegin;
      PetscSynchronizedPrintf(PETSC_COMM_WORLD,"Start Assembly\n");
      ierr = PetscTime(&v1); CHKERRQ(ierr);
      ierr = PetscGetCPUTime(&t1); CHKERRQ(ierr);
      PetscFunctionReturn(0);
    }
    
    PetscErrorCode postProcess() {
      PetscFunctionBegin;
      // Note MAT_FLUSH_ASSEMBLY
      ierr = MatAssemblyBegin(Aij,MAT_FLUSH_ASSEMBLY); CHKERRQ(ierr);
      ierr = MatAssemblyEnd(Aij,MAT_FLUSH_ASSEMBLY); CHKERRQ(ierr);
      ierr = VecAssemblyBegin(F); CHKERRQ(ierr);
      ierr = VecAssemblyEnd(F); CHKERRQ(ierr);
      ierr = PetscTime(&v2); CHKERRQ(ierr);
      ierr = PetscGetCPUTime(&t2); CHKERRQ(ierr);
      PetscSynchronizedPrintf(PETSC_COMM_WORLD,"End Assembly: Rank %d Time = %f CPU Time = %f\n",pcomm->rank(),v2-v1,t2-t1);
      PetscFunctionReturn(0);
    }
    
    virtual PetscErrorCode GetMatricesRows() {
      PetscFunctionBegin;
      //indicies ROWS
      //cout<<"GetMatricesRows " <<endl;
      row_mat = 0;
      ierr = GetRowGlobalIndices("CONCENTRATION",RowGlob[row_mat]); CHKERRQ(ierr);
      ierr = GetRowLocalIndices("CONCENTRATION",RowLocal[row_mat]); CHKERRQ(ierr);
      ierr = GetGaussRowNMatrix("CONCENTRATION",rowNMatrices[row_mat]); CHKERRQ(ierr);
      ierr = GetGaussRowDiffNMatrix("CONCENTRATION",rowDiffNMatrices[row_mat]); CHKERRQ(ierr);
      
      //      cout<<"RowLocal "<< (RowLocal[0])[0]<<"\t"<<(RowLocal[0])[1]<<"\t"<<(RowLocal[0])[2]<<"\t"<<(RowLocal[0])[3]<<"\t"<<endl;
      //      cout<<"RowGlob  "<< (RowGlob[0])[0]<<"\t"<<(RowGlob[0])[1]<<"\t"<<(RowGlob[0])[2]<<"\t"<<(RowGlob[0])[3]<<"\t"<<(RowGlob[0])[4]<<"\t"<<(RowGlob[0])[5]<<"\t"<<endl;
      
      row_mat++;
      for(int ee = 0;ee<6;ee++) { //edges matrices
        ierr = GetRowGlobalIndices("CONCENTRATION",MBEDGE,RowGlob[row_mat],ee); CHKERRQ(ierr);
        ierr = GetRowLocalIndices("CONCENTRATION",MBEDGE,RowLocal[row_mat],ee); CHKERRQ(ierr);
        if(RowGlob[row_mat].size()!=0) {
          ierr = GetGaussRowNMatrix("CONCENTRATION",MBEDGE,rowNMatrices[row_mat],ee); CHKERRQ(ierr);
          ierr = GetGaussRowDiffNMatrix("CONCENTRATION",MBEDGE,rowDiffNMatrices[row_mat],ee); CHKERRQ(ierr);
          row_mat++;
        }
      }
      for(int ff = 0;ff<4;ff++) { //faces matrices
        ierr = GetRowGlobalIndices("CONCENTRATION",MBTRI,RowGlob[row_mat],ff); CHKERRQ(ierr);
        ierr = GetRowLocalIndices("CONCENTRATION",MBTRI,RowLocal[row_mat],ff); CHKERRQ(ierr);
        if(RowGlob[row_mat].size()!=0) {
          ierr = GetGaussRowNMatrix("CONCENTRATION",MBTRI,rowNMatrices[row_mat],ff); CHKERRQ(ierr);
          ierr = GetGaussRowDiffNMatrix("CONCENTRATION",MBTRI,rowDiffNMatrices[row_mat],ff); CHKERRQ(ierr);
          row_mat++;
        }
      }
      ierr = GetRowGlobalIndices("CONCENTRATION",MBTET,RowGlob[row_mat]); CHKERRQ(ierr);
      ierr = GetRowLocalIndices("CONCENTRATION",MBTET,RowLocal[row_mat]); CHKERRQ(ierr);
      if(RowGlob[row_mat].size() != 0) { //volume matrices
        ierr = GetGaussRowNMatrix("CONCENTRATION",MBTET,rowNMatrices[row_mat]); CHKERRQ(ierr);
        ierr = GetGaussRowDiffNMatrix("CONCENTRATION",MBTET,rowDiffNMatrices[row_mat]); CHKERRQ(ierr);
        row_mat++;
      }
      PetscFunctionReturn(0);
    }
    
    virtual PetscErrorCode GetMatricesCols() {
      PetscFunctionBegin;
      //indicies COLS
      col_mat = 0;
      ierr = GetColGlobalIndices("CONCENTRATION",ColGlob[col_mat]); CHKERRQ(ierr);
      ierr = GetColLocalIndices("CONCENTRATION",ColLocal[col_mat]); CHKERRQ(ierr);
      ierr = GetGaussColNMatrix("CONCENTRATION",colNMatrices[col_mat]); CHKERRQ(ierr);
      ierr = GetGaussColDiffNMatrix("CONCENTRATION",colDiffNMatrices[col_mat]); CHKERRQ(ierr);
      //cout<<"colDiffNMatrices "<< (colDiffNMatrices[0])[0] <<endl;
      col_mat++;
      for(int ee = 0;ee<6;ee++) { //edges matrices
        ierr = GetColGlobalIndices("CONCENTRATION",MBEDGE,ColGlob[col_mat],ee); CHKERRQ(ierr);
        ierr = GetColLocalIndices("CONCENTRATION",MBEDGE,ColLocal[col_mat],ee); CHKERRQ(ierr);
        if(ColGlob[col_mat].size()!=0) {
          ierr = GetGaussColNMatrix("CONCENTRATION",MBEDGE,colNMatrices[col_mat],ee); CHKERRQ(ierr);
          ierr = GetGaussColDiffNMatrix("CONCENTRATION",MBEDGE,colDiffNMatrices[col_mat],ee); CHKERRQ(ierr);
          //cout<<"colDiffNMatrices Edges "<< (colDiffNMatrices[0])[0] <<endl;
          col_mat++;
        }
      }
      for(int ff = 0;ff<4;ff++) { //faces matrices
        ierr = GetColGlobalIndices("CONCENTRATION",MBTRI,ColGlob[col_mat],ff); CHKERRQ(ierr);
        ierr = GetColLocalIndices("CONCENTRATION",MBTRI,ColLocal[col_mat],ff); CHKERRQ(ierr);
        if(ColGlob[col_mat].size()!=0) {
          ierr = GetGaussColNMatrix("CONCENTRATION",MBTRI,colNMatrices[col_mat],ff); CHKERRQ(ierr);
          ierr = GetGaussColDiffNMatrix("CONCENTRATION",MBTRI,colDiffNMatrices[col_mat],ff); CHKERRQ(ierr);
          col_mat++;
        }
      }
      ierr = GetColGlobalIndices("CONCENTRATION",MBTET,ColGlob[col_mat]); CHKERRQ(ierr);
      ierr = GetColLocalIndices("CONCENTRATION",MBTET,ColLocal[col_mat]); CHKERRQ(ierr);
      if(ColGlob[col_mat].size() != 0) { //volume matrices
        ierr = GetGaussColNMatrix("CONCENTRATION",MBTET,colNMatrices[col_mat]); CHKERRQ(ierr);
        ierr = GetGaussColDiffNMatrix("CONCENTRATION",MBTET,colDiffNMatrices[col_mat]); CHKERRQ(ierr);
        col_mat++;
      }
      PetscFunctionReturn(0);
    }
    
    virtual PetscErrorCode GetMatrices() {
      PetscFunctionBegin;
      ierr = GetMatricesRows(); CHKERRQ(ierr);
      ierr = GetMatricesCols(); CHKERRQ(ierr);
      PetscFunctionReturn(0);
    }
    
    
    ublas::matrix<ublas::matrix<FieldData> > K;
    ublas::matrix<FieldData> BD;
    virtual PetscErrorCode Stiffness() {
      PetscFunctionBegin;
      double Moist_diffusivity;
      ierr = GetMatParameters(&Moist_diffusivity); CHKERRQ(ierr);
      //  cout<< "Moist_diffusivity "<<Moist_diffusivity<<endl;
      ierr = calculateD(Moist_diffusivity); CHKERRQ(ierr);
      //    cout<< "D "<<D<<endl;
      K.resize(row_mat,col_mat);
      int g_dim = g_NTET.size()/4;
      for(int rr = 0;rr<row_mat;rr++) {
        
        for(int gg = 0;gg<g_dim;gg++) {
          ublas::matrix<FieldData> &row_Mat = (rowDiffNMatrices[rr])[gg];
          double w = V*G_W_TET[gg];
          //        cout<<"(rowDiffNMatrices[rr])[gg]"<<(rowDiffNMatrices[rr])[gg]<<endl;
          BD.resize(3,row_Mat.size2());
          //ublas::noalias(BD) = prod( w*D,row_Mat );
          cblas_dsymm(CblasRowMajor,CblasLeft,CblasUpper,
                      BD.size1(),BD.size2(),
                      w,&*D.data().begin(),D.size2(),
                      &*row_Mat.data().begin(),row_Mat.size2(),
                      0.,&*BD.data().begin(),BD.size2());
          //cout<<"BD "<<BD<<endl;
          
          for(int cc = rr;cc<col_mat;cc++) {
            ublas::matrix<FieldData> &col_Mat = (rowDiffNMatrices[cc])[gg];
            //for first gauss point k= BT*D*B + 0*k
            // for gauss other than first  k= BT*D*B + 1*k
            double s = 1;
            if(gg == 0) {
              s = 0;
              K(rr,cc).resize(BD.size2(),col_Mat.size2());
            }
            //ublas::noalias(K(rr,cc)) = prod(trans(BD) , col_Mat ); // int BT*D*B
            cblas_dgemm(CblasRowMajor,CblasTrans,CblasNoTrans,
                        BD.size2(),col_Mat.size2(),BD.size1(),
                        1.,&*BD.data().begin(),BD.size2(),
                        &*col_Mat.data().begin(),col_Mat.size2(),
                        s,&*K(rr,cc).data().begin(),K(rr,cc).size2());
          }
        }
      }
      
      PetscFunctionReturn(0);
    }
    
    
    
    virtual PetscErrorCode Lhs() {
      PetscFunctionBegin;
      ierr = Stiffness(); CHKERRQ(ierr);
      for(int rr = 0;rr<row_mat;rr++) {
        if(RowGlob[rr].size()==0) continue;
        for(int cc = rr;cc<col_mat;cc++) {
          if(ColGlob[cc].size()==0) continue;
          if(RowGlob[rr].size()!=K(rr,cc).size1()) SETERRQ(PETSC_COMM_SELF,1,"data inconsistency");
          if(ColGlob[cc].size()!=K(rr,cc).size2()) SETERRQ(PETSC_COMM_SELF,1,"data inconsistency");
          ierr = MatSetValues(Aij,RowGlob[rr].size(),&(RowGlob[rr])[0],ColGlob[cc].size(),&(ColGlob[cc])[0],&(K(rr,cc).data())[0],ADD_VALUES); CHKERRQ(ierr);
          if(rr!=cc) {
            K(cc,rr) = trans(K(rr,cc));
            ierr = MatSetValues(Aij,ColGlob[cc].size(),&(ColGlob[cc])[0],RowGlob[rr].size(),&(RowGlob[rr])[0],&(K(cc,rr).data())[0],ADD_VALUES); CHKERRQ(ierr);
          }
          
        }
      }
      PetscFunctionReturn(0);
    }
    
    
    
    PetscErrorCode operator()() {
      PetscFunctionBegin;
      ierr = OpStudentStart_TET(g_NTET); CHKERRQ(ierr);
      
      //      cout<<"Hi from Moisture Transport Class "<<endl;
      ierr = GetMatrices(); CHKERRQ(ierr);
      ierr = Lhs(); CHKERRQ(ierr);
      
      ierr = OpStudentEnd(); CHKERRQ(ierr);
      PetscFunctionReturn(0);
    }
    
  };
  
  
}

#endif //__MOISTUREFEMETHOD_HPP__
